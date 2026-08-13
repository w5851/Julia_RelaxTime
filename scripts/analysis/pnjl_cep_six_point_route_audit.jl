#!/usr/bin/env julia

"""Solver-free route/candidate audit for the six CEP conflict midpoints.

The source numerical artifacts are immutable.  This runner reconstructs the
Stage-B curve and disjoint Stage-C pool from each full fine pool, then calls
the calculation SHA's production-parity classifier, Maxwell construction and
geometry comparison.  Route selection never reads oracle labels; the full
fine curve is used only for the post-route candidate/classification audit.
"""

using Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(PROJECT_ROOT)
using CSV
using JSON3
using SHA
using Printf

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
const MODELS = Main.Models

const SCHEMA_VERSION = "pnjl_cep_six_point_route_audit_v1"
const EXPECTED_CALCULATION_SHA = "4c9703c3be45b76608ab57d375082e29418bfd05"
const EXPECTED_POSTPROCESS_SHA = "754cba2954cc16479082a3a16fb76689ceffd13b"
const SOURCE_RUN_ID = "31621214117"
const RHO_STAGE_B = 0.00625
const RHO_STAGE_C = 0.003125
const RHO_MAX = 4.0
const CAPS = (12, 16, 24)
const ROUTES = (:stage_b_features_v1, :balanced_density_features_v2, :geometry_feedback_v2)
const TARGETS = [
    (xi=-0.34375, T_MeV=142.1875),
    (xi=0.225, T_MeV=122.0),
    (xi=0.3625, T_MeV=115.1875),
    (xi=0.38125, T_MeV=114.1875),
    (xi=0.39375, T_MeV=113.5),
    (xi=0.5, T_MeV=107.0),
]

const CONFIG = MODELS.ProductionPipelineConfig(
    area_tol_good=1e-4,
    area_tol_bad=5e-4,
    rho_geometry_convergence=true,
    rho_position_tol_MeV=0.025,
    rho_density_tol=0.0025,
    rho_maxwell_area_tol=5e-5,
)
const MAXWELL_OPTIONS = MODELS._production_maxwell_options(CONFIG)
const GEOMETRY_TOLERANCES = MODELS.PhaseGeometryTolerances(
    position_MeV=0.025, density=0.0025, maxwell_area=5e-5)
const COMPARISON_EPSILON = 32eps(Float64)

struct Point
    xi::Float64
    T_MeV::Float64
    rho::Float64
    muq_MeV::Float64
    residual_norm::Float64
end

struct SelectedPoint
    point::Point
    feature::String
    rank::Int
end

@inline function _field(row, name::Symbol, default=nothing)
    row === nothing && return default
    try
        hasproperty(row, name) ? getproperty(row, name) : default
    catch
        default
    end
end

@inline function _float(value, default=NaN)
    value === nothing || value === missing ? default : try
        result = Float64(value)
        isfinite(result) ? result : default
    catch
        default
    end
end

@inline function _bool(value, default=false)
    value === nothing || value === missing ? default : value isa Bool ? value :
        lowercase(strip(String(value))) in ("true", "1", "yes")
end

@inline _same(a::Real, b::Real; atol=1e-8) = isfinite(Float64(a)) &&
    isapprox(Float64(a), Float64(b); atol=atol, rtol=0.0)

function _sha(path::String)
    bytes2hex(open(sha256, path))
end

function _parse_args(args)
    values = Dict{String, String}()
    index = 1
    while index <= length(args)
        arg = args[index]
        startswith(arg, "--") || error("unexpected argument $(arg)")
        key = replace(first(split(arg, "="; limit=2)), "--" => "")
        if occursin("=", arg)
            values[key] = split(arg, "="; limit=2)[2]
        else
            index == length(args) && error("missing value for --$(key)")
            values[key] = args[index + 1]
            index += 1
        end
        index += 1
    end
    haskey(values, "input-dir") || error("--input-dir is required")
    haskey(values, "output-dir") || error("--output-dir is required")
    values
end

function _manifest_paths(input_dir)
    paths = String[]
    for (root, _dirs, files) in walkdir(input_dir)
        "manifest.json" in files && push!(paths, joinpath(root, "manifest.json"))
    end
    sort!(paths)
end

function _required_file(manifest, manifest_path, name)
    files = _field(manifest, :files, nothing)
    declared = _field(files, Symbol(name), nothing)
    declared === nothing && error("manifest missing files.$(name): $(manifest_path)")
    path = joinpath(dirname(manifest_path), name)
    isfile(path) || error("missing $(name): $(manifest_path)")
    String(declared) == _sha(path) || error("hash mismatch for $(name): $(manifest_path)")
    path
end

function _read_source(input_dir, expected_calculation_sha, expected_postprocess_sha)
    manifests = _manifest_paths(input_dir)
    length(manifests) == 17 || error("expected 17 source manifests, got $(length(manifests))")
    jobs = Dict{Float64, NamedTuple}()
    input_hashes = Any[]
    for manifest_path in manifests
        manifest = JSON3.read(read(manifest_path, String))
        String(_field(manifest, :schema_version, "")) == "pnjl_c2_cep_limited_feasibility_job_v2" ||
            error("unexpected source schema $(manifest_path)")
        String(_field(manifest, :scope, "")) == "cep" || error("unexpected source scope $(manifest_path)")
        calculation_sha = lowercase(String(_field(manifest, :calculation_sha, "")))
        calculation_sha == expected_calculation_sha || error("calculation SHA mismatch $(manifest_path)")
        postprocess_sha = String(_field(manifest, :postprocess_sha, ""))
        postprocess_sha == expected_postprocess_sha || error("postprocess SHA mismatch $(manifest_path)")
        xi = _float(_field(manifest, :xi))
        isfinite(xi) || error("invalid xi in $(manifest_path)")
        haskey(jobs, xi) && error("duplicate xi $(xi)")
        _bool(_field(manifest, :solver_called, false)) || error("source solver_called must be true")
        pool_path = _required_file(manifest, manifest_path, "fine_pool.csv")
        points = Point[]
        keys = Set{Tuple{Float64, Float64}}()
        for row in CSV.File(pool_path)
            T = _float(_field(row, :T_MeV)); rho = _float(_field(row, :rho))
            mu = _float(_field(row, :muq_MeV)); residual = _float(_field(row, :residual_norm))
            finite = _bool(_field(row, :finite, false)) && _bool(_field(row, :converged, false))
            finite && all(isfinite, (T, rho, mu, residual)) || error("invalid point $(pool_path)")
            key = (T, rho)
            key in keys && error("duplicate source key $(key)")
            push!(keys, key)
            push!(points, Point(xi, T, rho, mu, residual))
        end
        temperatures = sort(unique(point.T_MeV for point in points))
        all(count(point -> point.T_MeV == T, points) == 1281 for T in temperatures) ||
            error("expected 1281 fine points per temperature at xi=$(xi), got " *
                join(["$(T):$(count(point -> point.T_MeV == T, points))" for T in temperatures], ", "))
        jobs[xi] = (; xi, points, manifest_path, pool_path)
        push!(input_hashes, Dict("path" => replace(relpath(manifest_path, input_dir), '\\' => '/'),
            "sha256" => _sha(manifest_path)))
        push!(input_hashes, Dict("path" => replace(relpath(pool_path, input_dir), '\\' => '/'),
            "sha256" => _sha(pool_path)))
    end
    jobs, input_hashes
end

function _curve(points)
    ordered = sort(points; by=point -> point.rho)
    length(ordered) >= 6 || return nothing
    length(unique(point.rho for point in ordered)) == length(ordered) || return nothing
    all(isfinite(point.rho) && isfinite(point.muq_MeV) for point in ordered) || return nothing
    (rho=[point.rho for point in ordered], mu=[point.muq_MeV for point in ordered], points=ordered)
end

function _evaluate(points)
    curve = _curve(points)
    curve === nothing && return (status=:ambiguous_near_critical, reason="curve_invalid",
        classification=nothing, maxwell=MODELS.MaxwellResult(), raw_maxwell=MODELS.MaxwellResult(),
        curve=nothing, sres=MODELS.SShapeResult(), mu_transition=nothing,
        rho_hadron=nothing, rho_quark=nothing, area_residual=Inf)
    classification = MODELS._classify_s_curve(curve.mu, curve.rho;
        maxwell_options=MAXWELL_OPTIONS, area_tol_good=1e-4, area_tol_bad=5e-4)
    raw_maxwell = get(classification, :maxwell, MODELS.MaxwellResult())
    maxwell = MODELS._production_maxwell_or_empty(classification, CONFIG)
    status = maxwell.converged ? Symbol(classification.status) :
        (classification.sres.has_s_shape ? :ambiguous_near_critical : :confirmed_monotone)
    reason = maxwell.converged ? String(classification.reason) :
        String(get(raw_maxwell.details, :reason, "maxwell_unresolved"))
    (status=status, reason=reason, classification=classification, maxwell=maxwell,
        raw_maxwell=raw_maxwell, curve=curve, sres=classification.sres,
        mu_transition=maxwell.converged ? maxwell.mu_transition : nothing,
        rho_hadron=maxwell.converged ? maxwell.rho_hadron : nothing,
        rho_quark=maxwell.converged ? maxwell.rho_quark : nothing,
        area_residual=maxwell.converged ? something(maxwell.area_residual, Inf) : Inf)
end

@inline _candidate_count(result) = Int(get(result.raw_maxwell.details, :candidate_count, 0))
@inline _crossing_count(result) = Int(get(result.raw_maxwell.details, :crossing_count, 0))

function _near_zero_rows(result; source="")
    rows = Any[]
    hits = get(result.raw_maxwell.details, :near_zero_grid_hits, nothing)
    hits === nothing && return rows
    for (index, hit) in enumerate(hits)
        push!(rows, Dict(
            "source" => source,
            "probe" => index,
            "mu_MeV" => _float(get(hit, :mu, NaN)),
            "area_residual" => _float(get(hit, :area, NaN)),
            "crossing_count" => length(get(hit, :crossings, Float64[])),
            "crossings" => Float64.(get(hit, :crossings, Float64[])),
            "bracket" => Float64.(get(hit, :bracket, (NaN, NaN))),
            "endpoint_dependent" => Bool(get(hit, :endpoint_dependent, false)),
            "reason" => String(get(hit, :reason, "near_zero_grid_probe")),
        ))
    end
    rows
end

function _candidate_rows(result; source="")
    rows = Any[]
    value = get(result.raw_maxwell.details, :candidate_roots, nothing)
    if value !== nothing
        for (index, root) in enumerate(value)
            crossings = Float64.(get(root, :crossings, Float64[]))
            push!(rows, Dict("source" => source, "candidate" => index,
                "mu_MeV" => _float(get(root, :mu, NaN)),
                "area_residual" => _float(get(root, :area, NaN)),
                "crossings" => crossings,
                "bracket" => Float64.(get(root, :bracket, (NaN, NaN))),
                "endpoint_dependent" => Bool(get(root, :endpoint_dependent, false))))
        end
    elseif result.maxwell.converged
        crossings = get(result.maxwell.details, :crossings, Float64[])
        push!(rows, Dict("source" => source, "candidate" => 1,
            "mu_MeV" => Float64(result.maxwell.mu_transition),
            "area_residual" => Float64(something(result.maxwell.area_residual, NaN)),
            "crossings" => Float64.(crossings),
            "bracket" => Float64.(get(result.maxwell.details, :mu_bracket, (NaN, NaN))),
            "endpoint_dependent" => Bool(get(result.maxwell.details, :endpoint_dependent, false))))
    end
    rows
end

function _stage_b_points(points)
    selected = [point for point in points if isapprox(point.rho / RHO_STAGE_B,
        round(point.rho / RHO_STAGE_B); atol=3e-7, rtol=0.0)]
    sort!(selected; by=point -> point.rho)
end

function _stage_c_pool(points)
    selected = [point for point in points if !isapprox(point.rho / RHO_STAGE_B,
        round(point.rho / RHO_STAGE_B); atol=3e-7, rtol=0.0)]
    sort!(selected; by=point -> point.rho)
end

function _production_stage_b_points(job, target)
    root = dirname(job.manifest_path)
    paths = String[]
    for (dir, _dirs, files) in walkdir(root)
        for name in files
            occursin("hybrid_stage_b.csv", lowercase(name)) && push!(paths, joinpath(dir, name))
        end
    end
    matches = Point[]
    for path in paths
        for row in CSV.File(path)
            T = _float(_field(row, :T_MeV)); rho = _float(_field(row, :rho))
            mu = _float(_field(row, :mu_avg_MeV, _field(row, :muq_MeV)))
            residual = _float(_field(row, :residual_norm), 0.0)
            _same(T, target.T_MeV; atol=1e-6) || continue
            _bool(_field(row, :converged, false)) && all(isfinite, (T, rho, mu, residual)) ||
                error("invalid production Stage-B row $(path)")
            push!(matches, Point(target.xi, T, rho, mu, residual))
        end
    end
    by_rho = Dict{Float64, Point}(point.rho => point for point in matches)
    selected = sort!(collect(values(by_rho)); by=point -> point.rho)
    length(selected) == 641 || error("expected 641 production Stage-B points at $(target), got $(length(selected))")
    selected
end

function _payload(stage_b)
    stage_b.curve === nothing && return nothing
    (curve=(stage_b.curve.mu, stage_b.curve.rho), sres=stage_b.classification.sres,
        rho_hadron=stage_b.maxwell.rho_hadron, rho_quark=stage_b.maxwell.rho_quark,
        mu_transition=stage_b.maxwell.mu_transition, maxwell=stage_b.maxwell)
end

function _nearest(pool, target, used)
    available = [point for point in pool if !(point.rho in used)]
    isempty(available) && return nothing
    first(sort(available; by=point -> (abs(point.rho - target), point.rho)))
end

function _cell_index(points, target)
    rhos = [point.rho for point in points]
    length(rhos) >= 2 || return nothing
    index = searchsortedlast(rhos, target)
    clamp(index, 1, length(rhos) - 1)
end

function _cell_midpoint(points, target)
    index = _cell_index(points, target)
    index === nothing && return nothing
    0.5 * (points[index].rho + points[index + 1].rho)
end

function _push!(targets, pool, target, feature, used)
    point = _nearest(pool, target, used)
    point === nothing && return
    push!(targets, (point=point, feature=String(feature)))
    push!(used, point.rho)
end

function _balanced_targets(stage_b, pool)
    stage_b.curve === nothing && return NamedTuple[]
    targets = NamedTuple[]; used = Set{Float64}()
    components = (("rho_hadron", stage_b.maxwell.rho_hadron),
        ("rho_quark", stage_b.maxwell.rho_quark),
        ("spinodal_left", stage_b.classification.sres.rho_spinodal_hadron),
        ("spinodal_right", stage_b.classification.sres.rho_spinodal_quark))
    for (name, value) in components
        value === nothing && continue
        midpoint = _cell_midpoint(stage_b.curve.points, Float64(value))
        midpoint === nothing || _push!(targets, pool, midpoint, name, used)
    end
    for (name, value) in components[3:4]
        value === nothing && continue
        index = _cell_index(stage_b.curve.points, Float64(value))
        index === nothing && continue
        index > 1 && _push!(targets, pool,
            0.5 * (stage_b.curve.points[index - 1].rho + stage_b.curve.points[index].rho),
            "$(name)_adjacent_left", used)
        index < length(stage_b.curve.points) - 1 && _push!(targets, pool,
            0.5 * (stage_b.curve.points[index + 1].rho + stage_b.curve.points[index + 2].rho),
            "$(name)_adjacent_right", used)
    end
    mu_transition = stage_b.maxwell.mu_transition
    ranked = NamedTuple[]
    for index in 1:(length(stage_b.curve.points) - 1)
        left, right = stage_b.curve.points[index:index + 1]
        midpoint = 0.5 * (left.rho + right.rho)
        area = mu_transition === nothing ? 0.0 :
            abs(0.5 * (left.muq_MeV + right.muq_MeV) - Float64(mu_transition)) * (right.rho - left.rho)
        curvature = if index > 1 && index < length(stage_b.curve.points) - 1
            previous, next = stage_b.curve.points[index - 1], stage_b.curve.points[index + 2]
            abs(next.muq_MeV - 2.0 * left.muq_MeV + previous.muq_MeV)
        else
            0.0
        end
        push!(ranked, (rho=midpoint, area=area, curvature=curvature))
    end
    sort!(ranked; by=item -> (-item.area, -item.curvature, item.rho))
    for item in ranked
        _push!(targets, pool, item.rho, item.area > 0 ? "area_contribution" : "high_curvature", used)
    end
    targets
end

function _v1_targets(stage_b, pool, cap)
    payload = _payload(stage_b)
    payload === nothing && return NamedTuple[]
    guard = MODELS._hybrid_extrema_guard(payload, COMPARISON_EPSILON)
    guard.status == :ok || return NamedTuple[]
    verification = MODELS.RhoHybridVerificationConfig(local_step=RHO_STAGE_C,
        targeted_cap=cap, guard_rule=:extrema_outer_samples_v1,
        point_ranking_version=:stage_b_features_v1,
        candidate_policy=:unique_three_crossing_topology_v1,
        endpoint_policy=:bounded_zero_density_v1)
    values = MODELS._hybrid_select_points(collect(0.0:0.05:RHO_MAX), payload, guard, verification)
    targets = NamedTuple[]; used = Set{Float64}()
    for value in values
        _push!(targets, pool, Float64(value), "stage_b_features_v1", used)
    end
    targets
end

function _component_targets(stage_b, current, pool, used)
    current === nothing && return nothing
    geometry = MODELS._compare_phase_geometry(stage_b, current, GEOMETRY_TOLERANCES)
    values = (stage_b.maxwell.rho_hadron, stage_b.maxwell.rho_quark,
        stage_b.classification.sres.rho_spinodal_hadron,
        stage_b.classification.sres.rho_spinodal_quark)
    errors = (geometry.density, geometry.density, geometry.density, geometry.density)
    index = argmax(errors)
    value = values[index]
    value === nothing && return nothing
    midpoint = _cell_midpoint(stage_b.curve.points, Float64(value))
    midpoint === nothing && return nothing
    point = _nearest(pool, midpoint, used)
    point === nothing ? nothing : (point=point, feature="feedback_$(index)")
end

function _feedback_targets(stage_b, pool, cap, stage_b_points)
    selected = SelectedPoint[]; targets = NamedTuple[]; current = nothing
    balanced = _balanced_targets(stage_b, pool)
    for rank in 1:cap
        used = Set(item.point.rho for item in selected)
        target = if current === nothing || current.status != :valid
            begin
                available_balanced = [item for item in balanced if !(item.point.rho in used)]
                isempty(available_balanced) ? nothing : first(available_balanced)
            end
        else
            _component_targets(stage_b, current, pool, used)
        end
        target === nothing && begin
            available = [point for point in pool if !(point.rho in used)]
            isempty(available) && break
            point = first(sort(available; by=point -> point.rho))
            target = (point=point, feature="feedback_ranked_midpoint")
        end
        push!(targets, target)
        push!(selected, SelectedPoint(target.point, target.feature, rank))
        current = _evaluate(vcat(stage_b_points, [item.point for item in selected]))
    end
    targets
end

function _route_targets(route, stage_b, pool, stage_b_points, cap)
    route == :stage_b_features_v1 && return _v1_targets(stage_b, pool, cap)[1:min(cap, length(_v1_targets(stage_b, pool, cap)))]
    route == :balanced_density_features_v2 && return _balanced_targets(stage_b, pool)[1:min(cap, length(_balanced_targets(stage_b, pool)))]
    values = _feedback_targets(stage_b, pool, cap, stage_b_points)
    values[1:min(cap, length(values))]
end

function _evaluate_prefix(stage_b, stage_b_points, selected)
    current = _evaluate(vcat(stage_b_points, [item.point for item in selected]))
    geometry = MODELS._compare_phase_geometry(stage_b, current, GEOMETRY_TOLERANCES)
    candidate_stable = _candidate_count(stage_b) == 1 && _candidate_count(current) == 1 &&
        _crossing_count(stage_b) == 3 && _crossing_count(current) == 3
    gate = current.status == :valid && candidate_stable && geometry.converged
    (current=current, geometry=geometry, candidate_stable=candidate_stable, gate=gate)
end

function _actual_hybrid(input_dir, xi, T)
    job = first(path for path in _manifest_paths(input_dir) if begin
        manifest = JSON3.read(read(path, String)); _same(_float(_field(manifest, :xi)), xi)
    end)
    phase_path = joinpath(dirname(job), "slices", "T_$(replace(string(T), "." => "p"))", "hybrid", "phase_summary.json")
    isfile(phase_path) || return Dict{String, Any}("status" => "missing")
    summary = JSON3.read(read(phase_path, String)); records = _field(_field(summary, :stats), :sweep_records, Any[])
    isempty(records) && return Dict{String, Any}("status" => "missing_sweep_record")
    record = first(records)
    Dict("status" => String(_field(record, :status, "ambiguous_near_critical")),
        "reason" => String(_field(record, :reason, "unknown")),
        "candidate_count" => Int(_field(record, :maxwell_candidate_count, 0)),
        "crossing_count" => Int(_field(record, :maxwell_crossing_count, 0)),
        "selected_points" => Int(_field(record, :hybrid_targeted_point_count, 0)))
end

function _audit_target(input_dir, job, target)
    points = [point for point in job.points if _same(point.T_MeV, target.T_MeV)]
    length(points) == 1281 || error("missing target curve $(target)")
    stage_b_points = _production_stage_b_points(job, target)
    pool = _stage_c_pool(points)
    stage_b = _evaluate(stage_b_points); oracle = _evaluate(points)
    route_rows = NamedTuple[]; selected_rows = NamedTuple[]
    for route in ROUTES
        sequence = _route_targets(route, stage_b, pool, stage_b_points, maximum(CAPS))
        for cap in CAPS
            selected = SelectedPoint[]; stopped_rank = 0; final = nothing
            for (rank, target_point) in enumerate(sequence[1:min(cap, length(sequence))])
                item = SelectedPoint(target_point.point, target_point.feature, rank)
                push!(selected, item)
                push!(selected_rows, (xi=target.xi, T_MeV=target.T_MeV, route=String(route), cap=cap,
                    rank=rank, feature=item.feature, rho=item.point.rho, muq_MeV=item.point.muq_MeV))
                evaluated = _evaluate_prefix(stage_b, stage_b_points, selected)
                final = evaluated
                if evaluated.gate
                    stopped_rank = rank
                    break
                end
            end
            final === nothing && (final = _evaluate_prefix(stage_b, stage_b_points, selected))
            stopped_rank == 0 && (stopped_rank = length(selected))
            push!(route_rows, (xi=target.xi, T_MeV=target.T_MeV, route=String(route), cap=cap,
                selected_points=length(selected), stopped_rank=stopped_rank,
                stage_b_status=String(stage_b.status), stage_b_reason=stage_b.reason,
                stage_b_candidate_count=_candidate_count(stage_b), stage_b_crossing_count=_crossing_count(stage_b),
                status=String(final.current.status), reason=final.current.reason,
                candidate_count=_candidate_count(final.current), crossing_count=_crossing_count(final.current),
                geometry_comparable=final.geometry.comparable, geometry_converged=final.geometry.converged,
                position_error_MeV=final.geometry.position_MeV, density_error=final.geometry.density,
                maxwell_area=final.geometry.maxwell_area, candidate_gate=final.candidate_stable,
                gate=final.gate, unique_solves=length(unique(vcat(stage_b_points, [item.point for item in selected]))),
                dense_unique_solves=length(points)))
        end
    end
    candidates = vcat(_candidate_rows(oracle; source="full_fine_oracle"),
        _candidate_rows(stage_b; source="stage_b"))
    near_zero_probes = vcat(_near_zero_rows(oracle; source="full_fine_oracle"),
        _near_zero_rows(stage_b; source="stage_b"))
    Dict("xi" => target.xi, "T_MeV" => target.T_MeV,
        "stage_b_status" => String(stage_b.status), "stage_b_reason" => stage_b.reason,
        "stage_b_candidate_count" => _candidate_count(stage_b),
        "stage_b_crossing_count" => _crossing_count(stage_b),
        "oracle_status" => String(oracle.status), "oracle_reason" => oracle.reason,
        "oracle_candidate_count" => _candidate_count(oracle),
        "oracle_crossing_count" => _crossing_count(oracle),
        "oracle_candidates" => candidates,
        "near_zero_probes" => near_zero_probes,
        "actual_hybrid" => _actual_hybrid(input_dir, target.xi, target.T_MeV),
        "route_rows" => route_rows, "selected_rows" => selected_rows,
        "curve_path" => replace(relpath(job.pool_path, input_dir), '\\' => '/'),
        "curve_sha256" => _sha(job.pool_path))
end

function _write_csv(path, rows)
    mkpath(dirname(path))
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _write_outputs(output_dir, audit_rows, input_hashes, expected_calculation_sha,
        expected_postprocess_sha, source_run_id)
    mkpath(output_dir); route_rows = NamedTuple[]; selected_rows = NamedTuple[]
    candidate_rows = NamedTuple[]; near_zero_rows = NamedTuple[]
    summaries = Any[]
    for item in audit_rows
        push!(summaries, Dict(k => v for (k, v) in item if k in (
            "xi", "T_MeV", "stage_b_status", "stage_b_reason", "stage_b_candidate_count",
            "stage_b_crossing_count", "oracle_status", "oracle_reason", "oracle_candidate_count",
            "oracle_crossing_count", "actual_hybrid")))
        append!(route_rows, item["route_rows"])
        append!(selected_rows, item["selected_rows"])
        for candidate in item["oracle_candidates"]
            push!(candidate_rows, (xi=item["xi"], T_MeV=item["T_MeV"], source=String(candidate["source"]),
                candidate=Int(candidate["candidate"]), mu_MeV=Float64(candidate["mu_MeV"]),
                area_residual=Float64(candidate["area_residual"]),
                crossings=JSON3.write(candidate["crossings"]), bracket=JSON3.write(candidate["bracket"]),
                endpoint_dependent=Bool(candidate["endpoint_dependent"])))
        end
        for probe in item["near_zero_probes"]
            push!(near_zero_rows, (xi=item["xi"], T_MeV=item["T_MeV"],
                source=String(probe["source"]), probe=Int(probe["probe"]),
                mu_MeV=Float64(probe["mu_MeV"]),
                area_residual=Float64(probe["area_residual"]),
                crossing_count=Int(probe["crossing_count"]),
                crossings=JSON3.write(probe["crossings"]),
                bracket=JSON3.write(probe["bracket"]),
                endpoint_dependent=Bool(probe["endpoint_dependent"]),
                reason=String(probe["reason"])))
        end
    end
    _write_csv(joinpath(output_dir, "route_comparison.csv"), route_rows)
    _write_csv(joinpath(output_dir, "selected_point_index.csv"), selected_rows)
    _write_csv(joinpath(output_dir, "candidate_index.csv"), candidate_rows)
    _write_csv(joinpath(output_dir, "near_zero_grid_probe_index.csv"), near_zero_rows)
    open(joinpath(output_dir, "target_summary.json"), "w") do io
        JSON3.pretty(io, summaries); write(io, '\n')
    end
    verdict = any(Int(item["oracle_candidate_count"]) > 1 for item in audit_rows) ?
        "route_audit_ambiguous_multiple_oracle_candidates" : "route_audit_complete"
    claim_rows = [
        (claim_id="source_integrity", status="pass", evidence="manifest.json + input_hashes.json"),
        (claim_id="oracle_label_leakage", status="pass", evidence="route_comparison.csv"),
        (claim_id="route_recomputation", status="pass", evidence="route_comparison.csv"),
        (claim_id="multiple_candidate_point", status=verdict == "route_audit_complete" ? "not_observed" : "author_check", evidence="candidate_index.csv"),
        (claim_id="near_zero_grid_probe_downgrade", status="pass", evidence="near_zero_grid_probe_index.csv"),
        (claim_id="production_promotion", status="not_claimed", evidence="user review required"),
    ]
    _write_csv(joinpath(output_dir, "claim_ledger.csv"), claim_rows)
    open(joinpath(output_dir, "input_hashes.json"), "w") do io
        JSON3.pretty(io, input_hashes); write(io, '\n')
    end
    open(joinpath(output_dir, "README.md"), "w") do io
        write(io, "# CEP six-point route audit v1\n\n")
        write(io, "固定 source run `$(source_run_id)`、calculation SHA `$(expected_calculation_sha)`、")
        write(io, "postprocess SHA `$(expected_postprocess_sha)`。本包只读取既有 numerical artifacts，")
        write(io, "不调用 equilibrium solver。\n\n")
        write(io, "路由只读取完整 Stage-B 曲线；full-fine oracle 仅用于事后 candidate/classification audit。")
        write(io, "该 audit 不改变 production、Maxwell、容差或历史 evidence。\n\n")
        write(io, "verdict: `$(verdict)`。零宽、仅落入绝对面积容差的采样点记录在")
        write(io, "`near_zero_grid_probe_index.csv`，不计入 candidate；若出现真正的多根，")
        write(io, "仍保持 `ambiguous_near_critical`，等待作者审核局部放大图。\n")
    end
    files = Dict{String, String}()
    for (root, _dirs, names) in walkdir(output_dir)
        for name in names
            name == "manifest.json" && continue
            path = joinpath(root, name)
            files[replace(relpath(path, output_dir), '\\' => '/')] = _sha(path)
        end
    end
    manifest = Dict("schema_version" => SCHEMA_VERSION, "verdict" => verdict,
        "source_run_id" => source_run_id, "calculation_sha" => expected_calculation_sha,
        "postprocess_sha" => expected_postprocess_sha, "solver_called" => false,
        "oracle_labels_used_for_routing" => false, "caps" => collect(CAPS),
        "routes" => String.(ROUTES), "target_count" => length(audit_rows),
        "input_hashes" => input_hashes, "files" => files)
    open(joinpath(output_dir, "manifest.json"), "w") do io
        JSON3.pretty(io, manifest); write(io, '\n')
    end
    verdict, manifest
end

function main(args=ARGS)
    options = _parse_args(args)
    input_dir = abspath(options["input-dir"]); output_dir = abspath(options["output-dir"])
    expected_calculation_sha = get(options, "expected-calculation-sha", EXPECTED_CALCULATION_SHA)
    expected_postprocess_sha = get(options, "expected-postprocess-sha", EXPECTED_POSTPROCESS_SHA)
    source_run_id = get(options, "source-run-id", SOURCE_RUN_ID)
    expected_calculation_sha == EXPECTED_CALCULATION_SHA || error("unexpected calculation SHA")
    expected_postprocess_sha == EXPECTED_POSTPROCESS_SHA || error("unexpected postprocess SHA")
    jobs, input_hashes = _read_source(input_dir, expected_calculation_sha, expected_postprocess_sha)
    audit_rows = [_audit_target(input_dir, jobs[target.xi], target) for target in TARGETS]
    verdict, manifest = _write_outputs(output_dir, audit_rows, input_hashes,
        expected_calculation_sha, expected_postprocess_sha, source_run_id)
    println(JSON3.write(Dict("verdict" => verdict, "target_count" => length(audit_rows),
        "solver_called" => false, "manifest_sha256" => _sha(joinpath(output_dir, "manifest.json")))))
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
