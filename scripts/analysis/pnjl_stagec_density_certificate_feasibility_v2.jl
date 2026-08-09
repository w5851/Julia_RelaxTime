#!/usr/bin/env julia

"""Solver-free Stage-C density-certificate feasibility replay v2.

The source artifacts are produced by the immutable calculation checkout.  This
script only reads those curves and calls the phase classification/Maxwell/
geometry helpers; it never invokes the equilibrium solver or a production
pipeline.  The output is a versioned diagnostic package suitable for an
Actions aggregate replay.
"""

using Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(PROJECT_ROOT)

using CSV
using JSON3
using SHA
using Statistics
using Printf

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const MODELS = Main.Models

module StageCDensityCertificateFeasibilityV2

using CSV
using JSON3
using SHA
using Statistics
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MODELS = Main.Models
const SCHEMA_VERSION = "cep_hybrid_stagec_density_certificate_feasibility_v2"
const SOURCE_RUN_ID = "31296511813"
const SOURCE_CALCULATION_SHA = "ffa816df0a145f73d7490db1ed9ff10c92e017a4"
const SOURCE_POSTPROCESS_SHA = "4db46e1ef1694b28171ff79d6c00700a507b35ce"

const DENSITY_ANCHORS = [
    (-0.35, 51.0), (-0.25, 41.0), (-0.2, 41.0), (-0.15, 41.0),
    (-0.1, 41.0), (0.0, 51.0), (0.3, 21.0), (0.35, 51.0), (0.35, 101.0),
]
const FIRST_ORDER_CONTROLS = [(-0.5, 60.0), (0.0, 60.0), (0.5, 60.0)]
const MONOTONE_CONTROLS = [(-0.5, 160.0), (0.0, 145.0), (0.5, 120.0)]
const ALL_ANCHORS = vcat(DENSITY_ANCHORS, FIRST_ORDER_CONTROLS, MONOTONE_CONTROLS)
const METHODS = ("production_hybrid", "memoized_dense", "independent_oracle")
const ROUTES = (:stage_b_features_v1, :balanced_density_features_v2, :geometry_feedback_v2)
const CAPS = (12, 16, 24)

const STAGE_A_COARSE = 0.05
const STAGE_A_FINE = 0.025
const STAGE_B_COARSE = 0.0125
const STAGE_B_FINE = 0.00625
const STAGE_C_FINE = 0.003125
const POSITION_TOL = 0.025
const DENSITY_TOL = 0.0025
const AREA_TOL = 5e-5
const AREA_TOL_GOOD = 1e-4
const AREA_TOL_BAD = 5e-4
const FEATURE_RADIUS = 0.025
const COMPARISON_EPSILON = 32eps(Float64)

const PRODUCTION_CONFIG = MODELS.ProductionPipelineConfig(
    area_tol_good=AREA_TOL_GOOD,
    area_tol_bad=AREA_TOL_BAD,
    rho_geometry_convergence=true,
    rho_position_tol_MeV=POSITION_TOL,
    rho_density_tol=DENSITY_TOL,
    rho_maxwell_area_tol=AREA_TOL,
)
const PRODUCTION_MAXWELL_OPTIONS = MODELS._production_maxwell_options(PRODUCTION_CONFIG)
const GEOMETRY_TOLERANCES = MODELS.PhaseGeometryTolerances(
    position_MeV=POSITION_TOL,
    density=DENSITY_TOL,
    maxwell_area=AREA_TOL,
)

struct Point
    xi::Float64
    T_MeV::Float64
    rho::Float64
    muq_MeV::Float64
    residual_norm::Float64
    method::String
    level::Int
end

Point(xi::Real, T_MeV::Real, rho::Real, muq_MeV::Real,
    residual_norm::Real, method::AbstractString) = Point(
        Float64(xi), Float64(T_MeV), Float64(rho), Float64(muq_MeV),
        Float64(residual_norm), String(method), -1)

struct SelectedPoint
    point::Point
    feature::String
    rank::Int
    batch::Int
end

struct SourceData
    points::Vector{Point}
    slices::Vector{Any}
    costs::Vector{Any}
    manifests::Vector{Any}
    input_files::Vector{Any}
    job_count::Int
    source_run_conclusion::String
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
    value === nothing && return default
    value === missing && return default
    try
        result = Float64(value)
        isfinite(result) ? result : default
    catch
        default
    end
end

@inline function _bool(value, default=false)
    value === nothing && return default
    value === missing && return default
    value isa Bool && return value
    lowercase(strip(String(value))) in ("true", "1", "yes")
end

@inline function _json_field(value, name::Symbol, default=nothing)
    value === nothing && return default
    try
        hasproperty(value, name) ? getproperty(value, name) : default
    catch
        default
    end
end

@inline function _same(value::Real, target::Real; atol=2e-4)
    isfinite(Float64(value)) && isapprox(Float64(value), Float64(target); atol=atol, rtol=0.0)
end

@inline function _on_grid(rho::Float64, step::Float64)
    isapprox(rho / step, round(rho / step); atol=3e-7, rtol=0.0)
end

function _sha256_file(path::String)
    bytes2hex(open(sha256, path))
end

function _manifest_paths(input_dir::String)
    paths = String[]
    for (root, _dirs, files) in walkdir(input_dir)
        "manifest.json" in files || continue
        push!(paths, joinpath(root, "manifest.json"))
    end
    sort!(paths)
end

function _required_file(manifest, manifest_path::String, name::String)
    value = _json_field(_json_field(manifest, :files, nothing), Symbol(name), nothing)
    value === nothing && error("job manifest missing declared hash for $(name): $(manifest_path)")
    path = joinpath(dirname(manifest_path), name)
    isfile(path) || error("job artifact missing $(name): $(manifest_path)")
    actual = _sha256_file(path)
    String(value) == actual || error("job artifact hash mismatch for $(name): $(manifest_path)")
    path
end

function _read_job(manifest_path::String, expected_calculation_sha::String,
        expected_postprocess_sha::String, input_root::String)
    manifest = JSON3.read(read(manifest_path, String))
    scope = String(_json_field(manifest, :scope, ""))
    scope == "density_feasibility" || return nothing
    method = String(_json_field(manifest, :method, ""))
    method in METHODS || error("unexpected Stage-C method $(method)")
    xi = _float(_json_field(manifest, :xi, nothing))
    isfinite(xi) || error("invalid xi in $(manifest_path)")
    calculation_sha = String(_json_field(manifest, :calculation_sha, ""))
    calculation_sha == expected_calculation_sha || error("calculation SHA mismatch in $(manifest_path)")
    postprocess_sha = String(_json_field(manifest, :postprocess_sha,
        _json_field(manifest, :workflow_head_sha, "")))
    isempty(expected_postprocess_sha) || postprocess_sha == expected_postprocess_sha ||
        error("source postprocess SHA mismatch in $(manifest_path)")

    paths = Dict(name => _required_file(manifest, manifest_path, name)
        for name in ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv"))
    points = Point[]
    for row in CSV.File(paths["curve_points.csv"])
        row_method = String(_field(row, :method, method))
        row_method == method || error("curve method mismatch in $(manifest_path)")
        x = _float(_field(row, :xi, xi))
        t = _float(_field(row, :T_MeV, NaN))
        rho = _float(_field(row, :rho, NaN))
        mu = _float(_field(row, :muq_MeV, NaN))
        residual = _float(_field(row, :residual_norm, Inf), Inf)
        level = Int(round(_float(_field(row, :rho_level, -1.0), -1.0)))
        _bool(_field(row, :converged, false)) && _bool(_field(row, :finite, false)) ||
            error("non-finite or non-converged source curve point in $(manifest_path)")
        all(isfinite, (x, t, rho, mu, residual)) || error("invalid source curve point in $(manifest_path)")
        push!(points, Point(x, t, rho, mu, residual, method, level))
    end
    slices = Any[collect(CSV.File(paths["slice_metrics.csv"]))...]
    costs = Any[collect(CSV.File(paths["method_costs.csv"]))...]
    input_files = Any[(path=replace(relpath(path, input_root), '\\' => '/'),
        sha256=_sha256_file(path)) for path in values(paths)]
    push!(input_files, (path=replace(relpath(manifest_path, input_root), '\\' => '/'),
        sha256=_sha256_file(manifest_path)))
    return (; manifest, method, xi, points, slices, costs, input_files)
end

function _load_source(input_dir::String, expected_calculation_sha::String,
        expected_postprocess_sha::String, source_run_conclusion::String)
    jobs = Any[]
    for path in _manifest_paths(input_dir)
        job = _read_job(path, expected_calculation_sha, expected_postprocess_sha, input_dir)
        job === nothing || push!(jobs, job)
    end
    length(jobs) == 30 || error("expected 30 density-feasibility job artifacts, got $(length(jobs))")
    keys = Set{Tuple{String, Float64}}()
    points = Point[]
    slices = Any[]
    costs = Any[]
    manifests = Any[]
    input_files = Any[]
    for job in jobs
        key = (job.method, job.xi)
        key in keys && error("duplicate method/xi job $(key)")
        push!(keys, key)
        append!(points, job.points)
        append!(slices, job.slices)
        append!(costs, job.costs)
        push!(manifests, job.manifest)
        append!(input_files, job.input_files)
    end
    expected_keys = Set((method, xi) for method in METHODS for xi in
        (-0.5, -0.35, -0.25, -0.2, -0.15, -0.1, 0.0, 0.3, 0.35, 0.5))
    keys == expected_keys || error("source method/xi matrix does not match the 30-job contract")
    return SourceData(points, slices, costs, manifests, input_files, length(jobs), source_run_conclusion)
end

function _points(data::SourceData, method::String, anchor::Tuple{Float64, Float64}, predicate;
        levels=nothing)
    selected = Dict{Float64, Point}()
    for point in data.points
        point.method == method || continue
        _same(point.xi, anchor[1]; atol=1e-8) && _same(point.T_MeV, anchor[2]; atol=1e-8) || continue
        levels === nothing || point.level in levels || continue
        predicate(point.rho) || continue
        previous = get(selected, point.rho, nothing)
        previous === nothing || point.residual_norm < previous.residual_norm || continue
        selected[point.rho] = point
    end
    sort!(collect(values(selected)); by=point -> point.rho)
end

function _curve(points::Vector{Point})
    length(points) >= 6 || return nothing
    ordered = sort(unique(points); by=point -> point.rho)
    length(ordered) == length(points) || return nothing
    all(isfinite(point.rho) && isfinite(point.muq_MeV) for point in ordered) || return nothing
    (rho=[point.rho for point in ordered], mu=[point.muq_MeV for point in ordered], points=ordered)
end

function _empty_eval(reason::String)
    (status=:invalid, reason=reason, mu_transition=nothing, rho_hadron=nothing,
        rho_quark=nothing, area_residual=Inf, sres=MODELS.SShapeResult(),
        maxwell=MODELS.MaxwellResult(), raw_maxwell=MODELS.MaxwellResult(), curve=nothing)
end

function _evaluate(curve)
    curve === nothing && return _empty_eval("solver_or_curve_failure")
    cres = MODELS._classify_s_curve(curve.mu, curve.rho;
        maxwell_options=PRODUCTION_MAXWELL_OPTIONS,
        area_tol_good=AREA_TOL_GOOD, area_tol_bad=AREA_TOL_BAD)
    raw_maxwell = get(cres, :maxwell, MODELS.MaxwellResult())
    maxwell = MODELS._production_maxwell_or_empty(cres, PRODUCTION_CONFIG)
    status = maxwell.converged ? Symbol(cres.status) :
        (cres.sres.has_s_shape ? :unknown : Symbol(cres.status))
    reason = maxwell.converged ? String(cres.reason) :
        String(get(raw_maxwell.details, :reason, "maxwell_solver_tolerance_not_met"))
    return (
        status=status,
        reason=reason,
        mu_transition=maxwell.converged ? maxwell.mu_transition : nothing,
        rho_hadron=maxwell.converged ? maxwell.rho_hadron : nothing,
        rho_quark=maxwell.converged ? maxwell.rho_quark : nothing,
        area_residual=maxwell.converged ? something(maxwell.area_residual, Inf) : Inf,
        sres=cres.sres,
        maxwell=maxwell,
        raw_maxwell=raw_maxwell,
        curve=curve,
    )
end

function _geometry(left, right)
    MODELS._compare_phase_geometry(left, right, GEOMETRY_TOLERANCES)
end

function _semantic(left, right, geometry)
    if left.status == :invalid && right.status == :invalid &&
            left.reason == "no_s_shape" && right.reason == "no_s_shape" && geometry.converged
        return :confirmed_monotone
    elseif left.status == :valid && right.status == :valid && geometry.converged
        return :confirmed_first_order
    end
    return :ambiguous_near_critical
end

function _candidate_count(result)
    Int(get(result.raw_maxwell.details, :candidate_count, 0))
end

function _crossing_count(result)
    Int(get(result.raw_maxwell.details, :crossing_count, 0))
end

function _candidate_interval(result)
    value = get(result.raw_maxwell.details, :rho_interval, nothing)
    value === nothing && return nothing
    try
        low, high = Float64(value[1]), Float64(value[2])
        isfinite(low) && isfinite(high) && high > low ? (low, high) : nothing
    catch
        nothing
    end
end

function _candidate_stable(left, right)
    _candidate_count(left) == 1 && _candidate_count(right) == 1 || return false
    _crossing_count(left) == 3 && _crossing_count(right) == 3 || return false
    first_interval, second_interval = _candidate_interval(left), _candidate_interval(right)
    first_interval === nothing && return false
    second_interval === nothing && return false
    max(first_interval[1], second_interval[1]) <=
        min(first_interval[2], second_interval[2]) + STAGE_B_COARSE
end

function _component_geometry(left, right, geometry)
    names = ("rho_hadron", "rho_quark", "rho_spinodal_hadron", "rho_spinodal_quark")
    if geometry.comparable && left.status == :invalid && right.status == :invalid &&
            left.reason == "no_s_shape" && right.reason == "no_s_shape"
        return [(component=name, error=0.0, normalized=0.0, comparable=true, pass=true)
            for name in names]
    end
    if geometry.comparable && left.status == :valid && right.status == :valid
        left_snapshot = MODELS._phase_geometry_snapshot(left)
        right_snapshot = MODELS._phase_geometry_snapshot(right)
        errors = abs.(Float64.(left_snapshot.density) .- Float64.(right_snapshot.density))
        return [(component=names[index], error=errors[index], normalized=errors[index] / DENSITY_TOL,
            comparable=true, pass=errors[index] <= DENSITY_TOL) for index in eachindex(names)]
    end
    [(component=name, error=Inf, normalized=Inf, comparable=false, pass=false) for name in names]
end

function _base_curve_payload(result)
    result.curve === nothing && return nothing
    (curve=(result.curve.mu, result.curve.rho), sres=result.sres,
        rho_hadron=result.rho_hadron, rho_quark=result.rho_quark,
        mu_transition=result.mu_transition, maxwell=result.maxwell)
end

function _pool_nearest(pool::Vector{Point}, target::Float64, used::Set{Float64})
    available = [point for point in pool if !(point.rho in used)]
    isempty(available) && return nothing
    first(sort!(available; by=point -> (abs(point.rho - target), point.rho)))
end

function _cell_index(points::Vector{Point}, target::Float64)
    length(points) >= 2 || return nothing
    rhos = [point.rho for point in points]
    index = searchsortedlast(rhos, target)
    clamp(index, 1, length(points) - 1)
end

function _cell_midpoint(points::Vector{Point}, target::Float64)
    index = _cell_index(points, target)
    index === nothing && return nothing
    0.5 * (points[index].rho + points[index + 1].rho)
end

function _push_target!(targets::Vector{Tuple{Float64, String}}, pool, target::Float64,
        feature::String, used::Set{Float64})
    point = _pool_nearest(pool, target, used)
    point === nothing && return
    push!(targets, (point.rho, feature))
    push!(used, point.rho)
end

function _balanced_targets(stage_b, pool::Vector{Point})
    stage_b.curve === nothing && return Tuple{Float64, String}[]
    targets = Tuple{Float64, String}[]
    used = Set{Float64}()
    components = (
        ("rho_hadron", stage_b.rho_hadron),
        ("rho_quark", stage_b.rho_quark),
        ("spinodal_left", stage_b.sres.rho_spinodal_hadron),
        ("spinodal_right", stage_b.sres.rho_spinodal_quark),
    )
    for (name, value) in components
        value === nothing && continue
        midpoint = _cell_midpoint(stage_b.curve.points, Float64(value))
        midpoint === nothing || _push_target!(targets, pool, midpoint, name, used)
    end
    for (name, value) in components[3:4]
        value === nothing && continue
        index = _cell_index(stage_b.curve.points, Float64(value))
        index === nothing && continue
        index > 1 && _push_target!(targets, pool,
            0.5 * (stage_b.curve.points[index - 1].rho + stage_b.curve.points[index].rho),
            "$(name)_adjacent_left", used)
        index < length(stage_b.curve.points) - 1 && _push_target!(targets, pool,
            0.5 * (stage_b.curve.points[index + 1].rho + stage_b.curve.points[index + 2].rho),
            "$(name)_adjacent_right", used)
    end

    mu_transition = stage_b.mu_transition
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
        _push_target!(targets, pool, item.rho, item.area > 0 ? "area_contribution" : "high_curvature", used)
    end
    targets
end

function _v1_targets(stage_b, pool::Vector{Point}, cap::Int)
    payload = _base_curve_payload(stage_b)
    payload === nothing && return Tuple{Float64, String}[]
    guard = MODELS._hybrid_extrema_guard(payload, COMPARISON_EPSILON)
    guard.status == :ok || return Tuple{Float64, String}[]
    verification = MODELS.RhoHybridVerificationConfig(
        local_step=STAGE_C_FINE, targeted_cap=cap,
        guard_rule=:extrema_outer_samples_v1,
        point_ranking_version=:stage_b_features_v1,
        candidate_policy=:unique_three_crossing_topology_v1,
        endpoint_policy=:bounded_zero_density_v1,
    )
    base_grid = collect(0.0:STAGE_A_COARSE:4.0)
    rho_values = MODELS._hybrid_select_points(base_grid, payload, guard, verification)
    targets = Tuple{Float64, String}[]
    used = Set{Float64}()
    for rho in rho_values
        point = _pool_nearest(pool, Float64(rho), used)
        point === nothing && continue
        push!(targets, (point.rho, "stage_b_features_v1"))
        push!(used, point.rho)
    end
    targets
end

function _feedback_target(stage_b, pool, selected::Vector{SelectedPoint}, current)
    used = Set(point.point.rho for point in selected)
    if current === nothing || current.status != :valid
        balanced = _balanced_targets(stage_b, pool)
        for target in balanced
            target[1] in used || return target
        end
    else
        geometry = _geometry(stage_b, current)
        components = _component_geometry(stage_b, current, geometry)
        worst = findmax([row.normalized for row in components])[2]
        target_value = (
            stage_b.rho_hadron, stage_b.rho_quark,
            stage_b.sres.rho_spinodal_hadron, stage_b.sres.rho_spinodal_quark,
        )[worst]
        if target_value !== nothing
            midpoint = _cell_midpoint(stage_b.curve.points, Float64(target_value))
            midpoint === nothing || begin
                point = _pool_nearest(pool, midpoint, used)
                point === nothing || return (point.rho, "feedback_$(components[worst].component)")
            end
        end
    end
    candidates = [point for point in pool if !(point.rho in used)]
    isempty(candidates) && return nothing
    point = first(sort(candidates; by=item -> item.rho))
    (point.rho, "feedback_ranked_midpoint")
end

function _next_target(route::Symbol, stage_b, pool, selected, current, cap)
    used = Set(item.point.rho for item in selected)
    targets = route == :stage_b_features_v1 ? _v1_targets(stage_b, pool, cap) :
        _balanced_targets(stage_b, pool)
    if route == :geometry_feedback_v2
        return _feedback_target(stage_b, pool, selected, current)
    end
    for target in targets
        target[1] in used || return target
    end
    nothing
end

function _candidate_targets(route::Symbol, stage_b, pool, stage_b_points, cap::Int)
    if route == :stage_b_features_v1
        targets = _v1_targets(stage_b, pool, cap)
        return targets[1:min(cap, length(targets))]
    end
    route == :balanced_density_features_v2 && begin
        targets = _balanced_targets(stage_b, pool)
        return targets[1:min(cap, length(targets))]
    end
    selected = SelectedPoint[]
    current = nothing
    targets = Tuple{Float64, String}[]
    for batch in 1:cap
        target = _feedback_target(stage_b, pool, selected, current)
        target === nothing && break
        rho, feature = target
        index = findfirst(point -> isapprox(point.rho, rho; atol=3e-7, rtol=0.0), pool)
        index === nothing && break
        point = pool[index]
        push!(targets, (point.rho, feature))
        push!(selected, SelectedPoint(point, feature, length(selected) + 1, batch))
        current = _evaluate(_curve(_merge_points(stage_b_points, selected)))
    end
    targets
end

function _merge_points(base::Vector{Point}, selected::Vector{SelectedPoint})
    by_rho = Dict(point.rho => point for point in base)
    for item in selected
        by_rho[item.point.rho] = item.point
    end
    sort!(collect(values(by_rho)); by=point -> point.rho)
end

function _certificate(stage_b_coarse, stage_b_fine, current)
    stage_b_geometry = _geometry(stage_b_coarse, stage_b_fine)
    stage_b_semantic = _semantic(stage_b_coarse, stage_b_fine, stage_b_geometry)
    if stage_b_semantic == :confirmed_monotone
        return (status=:confirmed_monotone, reason="stage_b_monotone_certificate",
            geometry=stage_b_geometry, stage_b_geometry=stage_b_geometry,
            component_geometry=_component_geometry(stage_b_fine, stage_b_fine, stage_b_geometry))
    end
    current === nothing && return (status=:ambiguous_near_critical, reason="stage_c_not_run",
        geometry=stage_b_geometry, stage_b_geometry=stage_b_geometry,
        component_geometry=_component_geometry(stage_b_fine, stage_b_fine, stage_b_geometry))
    current_geometry = _geometry(stage_b_fine, current)
    stable = _candidate_stable(stage_b_fine, current)
    valid = stage_b_fine.status == :valid && current.status == :valid &&
        stable && current_geometry.converged
    return (status=valid ? :confirmed_first_order : :ambiguous_near_critical,
        reason=valid ? "stage_c_recomputed_geometry_certificate" :
            (current.status == :unknown ? "stage_c_maxwell_unresolved" : "stage_c_geometry_unresolved"),
        geometry=current_geometry, stage_b_geometry=stage_b_geometry,
        component_geometry=_component_geometry(stage_b_fine, current, current_geometry))
end

function _oracle_status(data::SourceData, anchor)
    rows = [row for row in data.slices if String(_field(row, :method, "")) == "independent_oracle" &&
        _same(_float(_field(row, :xi, NaN)), anchor[1]; atol=1e-8) &&
        _same(_float(_field(row, :T_MeV, NaN)), anchor[2]; atol=1e-8)]
    isempty(rows) ? "missing" : String(_field(first(rows), :result_status, "missing"))
end

function _expected_status(anchor, oracle_status::String)
    anchor in FIRST_ORDER_CONTROLS && return "confirmed_first_order"
    anchor in MONOTONE_CONTROLS && return "confirmed_monotone"
    oracle_status
end

function _base_points(data::SourceData, anchor, stage_a_monotone::Bool,
        stage_a_fine::Vector{Point}, stage_b_fine::Vector{Point}=Point[])
    source = stage_a_monotone ? stage_a_fine : vcat(stage_a_fine, stage_b_fine)
    unique((point.xi, point.T_MeV, point.rho) for point in source)
end

function _run_anchor(data::SourceData, anchor, route::Symbol, cap::Int)
    # The numerical matrix deliberately separates the production Stage-A
    # route from the independent-oracle Stage-B/Stage-C session.  The
    # memoized-dense method is reserved for the cost baseline below.
    # rho_level is authoritative here; a numeric grid predicate alone would
    # accidentally admit targeted points that happen to lie on a base grid.
    stage_a_coarse_points = _points(data, "production_hybrid", anchor,
        rho -> _on_grid(rho, STAGE_A_COARSE); levels=(0,))
    stage_a_fine_points = _points(data, "production_hybrid", anchor,
        rho -> _on_grid(rho, STAGE_A_FINE); levels=(0, 1))
    stage_b_coarse_points = _points(data, "independent_oracle", anchor,
        rho -> _on_grid(rho, STAGE_B_COARSE); levels=(0,))
    stage_b_fine_points = _points(data, "independent_oracle", anchor,
        rho -> _on_grid(rho, STAGE_B_FINE); levels=(0,))
    stage_c_pool = _points(data, "independent_oracle", anchor,
        rho -> _on_grid(rho, STAGE_C_FINE) && !_on_grid(rho, STAGE_B_FINE); levels=(1,))
    isempty(intersect(Set(point.rho for point in stage_b_fine_points),
        Set(point.rho for point in stage_c_pool))) ||
        error("Stage-C pool overlaps the complete Stage-B curve at $(anchor)")
    stage_a_coarse, stage_a_fine = _evaluate(_curve(stage_a_coarse_points)), _evaluate(_curve(stage_a_fine_points))
    stage_a_geometry = _geometry(stage_a_coarse, stage_a_fine)
    stage_a_semantic = _semantic(stage_a_coarse, stage_a_fine, stage_a_geometry)
    if stage_a_semantic == :confirmed_monotone
        return (xi=anchor[1], T_MeV=anchor[2], route=route, cap=cap,
            oracle_status=_oracle_status(data, anchor), expected_status=_expected_status(anchor, _oracle_status(data, anchor)),
            result_status=:confirmed_monotone, reason="stage_a_monotone_certificate",
            stage_a_status=stage_a_semantic, stage_b_status=:not_run, stage_c_status=:not_run,
            stage_a_geometry=stage_a_geometry, stage_b_geometry=nothing, geometry=stage_a_geometry,
            components=_component_geometry(stage_a_fine, stage_a_fine, stage_a_geometry),
            selected=SelectedPoint[], base_keys=_base_points(data, anchor, true, stage_a_fine_points),
            candidate_targets=Tuple{Float64, String}[], candidate_count=0, crossing_count=0, finite=true)
    end

    stage_b_coarse, stage_b_fine = _evaluate(_curve(stage_b_coarse_points)), _evaluate(_curve(stage_b_fine_points))
    stage_b_payload = stage_b_fine
    candidate_targets = _candidate_targets(route, stage_b_payload, stage_c_pool,
        stage_b_fine_points, cap)
    stage_b_initial = _certificate(stage_b_coarse, stage_b_fine, nothing)
    if stage_b_initial.status == :confirmed_monotone ||
            (stage_b_fine.status == :valid && stage_b_initial.geometry.converged &&
             _candidate_stable(stage_b_coarse, stage_b_fine))
        status = stage_b_initial.status == :confirmed_monotone ? :confirmed_monotone : :confirmed_first_order
        return (xi=anchor[1], T_MeV=anchor[2], route=route, cap=cap,
            oracle_status=_oracle_status(data, anchor), expected_status=_expected_status(anchor, _oracle_status(data, anchor)),
            result_status=status, reason="stage_b_certificate", stage_a_status=stage_a_semantic,
            stage_b_status=status, stage_c_status=:not_run, stage_a_geometry=stage_a_geometry,
            stage_b_geometry=stage_b_initial.geometry, geometry=stage_b_initial.geometry,
            components=stage_b_initial.component_geometry, selected=SelectedPoint[],
            candidate_targets=candidate_targets,
            base_keys=_base_points(data, anchor, false, stage_a_fine_points, stage_b_fine_points),
            candidate_count=max(_candidate_count(stage_b_coarse), _candidate_count(stage_b_fine)),
            crossing_count=max(_crossing_count(stage_b_coarse), _crossing_count(stage_b_fine)), finite=true)
    end

    selected = SelectedPoint[]
    current = nothing
    for batch in 1:cap
        target = _next_target(route, stage_b_payload, stage_c_pool, selected, current, cap)
        target === nothing && break
        rho, feature = target
        point = findfirst(item -> isapprox(item.rho, rho; atol=3e-7, rtol=0.0), stage_c_pool)
        point === nothing && break
        push!(selected, SelectedPoint(stage_c_pool[point], feature, length(selected) + 1, batch))
        current = _evaluate(_curve(_merge_points(stage_b_fine_points, selected)))
        certificate = _certificate(stage_b_coarse, stage_b_fine, current)
        certificate.status == :confirmed_first_order && break
    end
    certificate = _certificate(stage_b_coarse, stage_b_fine, current)
    base_keys = _base_points(data, anchor, false, stage_a_fine_points, stage_b_fine_points)
    selected_keys = [(item.point.xi, item.point.T_MeV, item.point.rho) for item in selected]
    result = current === nothing ? stage_b_fine : current
    return (xi=anchor[1], T_MeV=anchor[2], route=route, cap=cap,
        oracle_status=_oracle_status(data, anchor), expected_status=_expected_status(anchor, _oracle_status(data, anchor)),
        result_status=certificate.status, reason=certificate.reason,
        stage_a_status=stage_a_semantic, stage_b_status=_semantic(stage_b_coarse, stage_b_fine, stage_b_initial.geometry),
        stage_c_status=certificate.status, stage_a_geometry=stage_a_geometry,
        stage_b_geometry=stage_b_initial.geometry, geometry=certificate.geometry,
        components=certificate.component_geometry, selected=selected,
        candidate_targets=candidate_targets,
        base_keys=vcat(base_keys, selected_keys), candidate_count=max(_candidate_count(stage_b_fine), _candidate_count(result)),
        crossing_count=max(_crossing_count(stage_b_fine), _crossing_count(result)), finite=true)
end

function _dense_cost(data::SourceData)
    sum(_float(_field(row, :unique_solves, 0.0), 0.0) for row in data.costs if
        String(_field(row, :method, "")) == "memoized_dense") |> round |> Int
end

function _frontier(data::SourceData)
    dense = _dense_cost(data)
    frontiers = NamedTuple[]
    all_rows = NamedTuple[]
    component_rows = NamedTuple[]
    selected_rows = NamedTuple[]
    candidate_rows = NamedTuple[]
    for route in ROUTES, cap in CAPS
        results = [_run_anchor(data, anchor, route, cap) for anchor in ALL_ANCHORS]
        mismatch = count(result -> result.result_status != Symbol(result.expected_status), results)
        unsupported = count(result -> result.oracle_status == "ambiguous_near_critical" &&
            result.result_status == :confirmed_first_order, results)
        multiple = count(result -> result.candidate_count > 1 || result.crossing_count > 3, results)
        geometry_failures = count(result -> result.result_status == :confirmed_first_order &&
            (result.geometry === nothing || !result.geometry.converged), results)
        finite_failures = count(result -> !result.finite, results)
        unique_solves = sum(length(unique(result.base_keys)) for result in results)
        targeted_points = sum(length(result.selected) for result in results)
        state_gate = mismatch == 0 && unsupported == 0
        geometry_gate = geometry_failures == 0
        candidate_gate = multiple == 0
        finite_gate = finite_failures == 0
        cost_gate = unique_solves <= dense
        feasible = state_gate && geometry_gate && candidate_gate && finite_gate && cost_gate
        push!(frontiers, (route=String(route), cap=cap, anchors=length(results),
            classification_mismatches=mismatch, unsupported_confirmations=unsupported,
            geometry_failures=geometry_failures, multiple_candidate_anchors=multiple,
            finite_failures=finite_failures, unique_solves=unique_solves,
            targeted_points=targeted_points, dense_unique_solves=dense,
            state_gate=state_gate, geometry_gate=geometry_gate, candidate_gate=candidate_gate,
            finite_gate=finite_gate, cost_gate=cost_gate, feasible=feasible))
        for result in results
            push!(all_rows, (xi=result.xi, T_MeV=result.T_MeV, route=String(route), cap=cap,
                oracle_status=result.oracle_status, expected_status=result.expected_status,
                simulated_status=String(result.result_status), reason=result.reason,
                stage_a_status=String(result.stage_a_status), stage_b_status=String(result.stage_b_status),
                stage_c_status=String(result.stage_c_status), candidate_count=result.candidate_count,
                crossing_count=result.crossing_count, selected_points=length(result.selected),
                geometry_evaluable=result.geometry === nothing ? false : result.geometry.comparable,
                geometry_pass=result.geometry === nothing ? false : result.geometry.converged,
                geometry_unresolved=result.geometry === nothing ? true : !result.geometry.comparable,
                position_error_MeV=result.geometry === nothing ? Inf : result.geometry.position_MeV,
                density_error=result.geometry === nothing ? Inf : result.geometry.density,
                maxwell_area_gate=result.geometry === nothing ? Inf : result.geometry.maxwell_area,
                unique_solves=length(unique(result.base_keys)), dense_unique_solves=dense, cost_gate=length(unique(result.base_keys)) <= dense,
                finite_gate=result.finite))
            for (index, item) in enumerate(result.selected)
                push!(selected_rows, (xi=result.xi, T_MeV=result.T_MeV, route=String(route), cap=cap,
                    rank=item.rank, batch=item.batch, feature=item.feature, rho=item.point.rho, muq_MeV=item.point.muq_MeV))
            end
            for (index, item) in enumerate(result.candidate_targets)
                push!(candidate_rows, (xi=result.xi, T_MeV=result.T_MeV, route=String(route), cap=cap,
                    rank=index, feature=item[2], rho=item[1]))
            end
            for row in result.components
                push!(component_rows, (xi=result.xi, T_MeV=result.T_MeV, route=String(route), cap=cap,
                    component=row.component, error=row.error, normalized_error=row.normalized,
                    comparable=row.comparable, pass=row.pass))
            end
        end
    end
    return (; dense, frontiers, all_rows, component_rows, selected_rows, candidate_rows)
end

function _verdict(frontier)
    cap12 = [row for row in frontier.frontiers if row.cap == 12 && row.feasible]
    cap_high = [row for row in frontier.frontiers if row.cap in (16, 24) && row.feasible]
    any(row -> row.finite_failures > 0, frontier.frontiers) && return "artifact_invalid"
    any(row -> row.unsupported_confirmations > 0, frontier.frontiers) && return "oracle_inconclusive"
    any(row -> row.multiple_candidate_anchors > 0, frontier.frontiers) && return "maxwell_candidate_inconclusive"
    !isempty(cap12) && return "feasible_candidate"
    !isempty(cap_high) && return "cap_contract_inconclusive"
    all(row -> row.cost_gate == false, [row for row in frontier.frontiers if row.cap == 12]) && return "performance_inconclusive"
    "integration_failed"
end

function _selected_policy(frontier, verdict)
    candidates = [row for row in frontier.frontiers if row.cap == 12 && row.feasible]
    priority = Dict("stage_b_features_v1" => 1, "balanced_density_features_v2" => 2, "geometry_feedback_v2" => 3)
    selected = isempty(candidates) ? nothing : first(sort(candidates; by=row ->
        (row.unique_solves, row.targeted_points, priority[row.route])))
    Dict{String, Any}(
        "schema_version" => SCHEMA_VERSION,
        "verdict" => verdict,
        "route" => selected === nothing ? nothing : selected.route,
        "cap" => selected === nothing ? nothing : selected.cap,
        "point_ranking_version" => selected === nothing ? nothing : selected.route,
        "stage_a_coarse_step" => STAGE_A_COARSE,
        "stage_b_fine_step" => STAGE_B_FINE,
        "stage_c_local_step" => STAGE_C_FINE,
        "target_cap" => 12,
        "feature_radius" => FEATURE_RADIUS,
        "position_tol_MeV" => POSITION_TOL,
        "density_tol" => DENSITY_TOL,
        "area_tol" => AREA_TOL,
        "dense_unique_solves" => frontier.dense,
        "reason" => selected === nothing ? "no cap-12 route satisfies all gates" : "minimum-cost cap-12 route",
    )
end

function _write_csv(path::String, rows)
    mkpath(dirname(path))
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _write_outputs(output_dir::String, data::SourceData, frontier, verdict, policy,
        expected_calculation_sha::String, expected_postprocess_sha::String, producer_sha::String)
    mkpath(output_dir)
    _write_csv(joinpath(output_dir, "route_comparison.csv"), frontier.all_rows)
    _write_csv(joinpath(output_dir, "component_geometry.csv"), frontier.component_rows)
    _write_csv(joinpath(output_dir, "selected_point_index.csv"), frontier.selected_rows)
    _write_csv(joinpath(output_dir, "candidate_point_index.csv"), frontier.candidate_rows)
    _write_csv(joinpath(output_dir, "cost_frontier.csv"), frontier.frontiers)
    (output_dir * "/selected_policy.json") |> path -> open(path, "w") do io
        JSON3.pretty(io, policy); write(io, '\n')
    end
    _write_csv(joinpath(output_dir, "author_adjudications.csv"), [
        (xi=xi, T_MeV=T, expected_status="confirmed_first_order", source="oracle_non_ambiguous_fixture")
        for (xi, T) in FIRST_ORDER_CONTROLS
    ])
    _write_csv(joinpath(output_dir, "claim_ledger.csv"), [
        (claim_id="classification", status=verdict == "feasible_candidate" ? "pass" : "inconclusive", evidence="route_comparison.csv"),
        (claim_id="geometry", status=all(row -> row.geometry_gate, frontier.frontiers) ? "pass" : "inconclusive", evidence="component_geometry.csv"),
        (claim_id="cost", status=all(row -> row.cost_gate, frontier.frontiers) ? "pass" : "inconclusive", evidence="cost_frontier.csv"),
        (claim_id="promotion", status="not_claimed", evidence="requires targeted/full shadow and author review"),
    ])
    open(joinpath(output_dir, "README.md"), "w") do io
        write(io, "# Stage-C density certificate feasibility v2\n\n")
        write(io, "verdict: `$(verdict)`。固定 source run `$(SOURCE_RUN_ID)` 的 solver-free replay；\n")
        write(io, "不调用 equilibrium solver，不修改 v1 evidence、production、reference 或 transport。\n\n")
        write(io, "三条 route 在合并 Stage-B 全域曲线和实际 Stage-C 点后重新执行 Julia Maxwell/geometry；\n")
        write(io, "Stage A 使用 production_hybrid；Stage B 与 Stage-C pool 均来自\n")
        write(io, "independent_oracle session 的 level=0/1；rho_level 是网格层级的权威字段。\n")
        write(io, "oracle 标签只用于事后 gate。cap 16/24 仅作诊断，不能授权 production。\n")
    end
    open(joinpath(output_dir, "AUDIT.md"), "w") do io
        write(io, "# Stage-C density certificate feasibility v2 audit\n\n")
        write(io, "成本按 Stage-A、完整 Stage-B 与实际 Stage-C rho key 并集计费；旧 production\n")
        write(io, "diagnostics 不参与 v2 geometry 证书。`solver_called=false`。\n")
    end
    open(joinpath(output_dir, "plot_manifest.json"), "w") do io
        JSON3.pretty(io, Dict("schema_version" => SCHEMA_VERSION,
            "figures" => Any[], "raw_curves_external" => true)); write(io, '\n')
    end
    # Write the manifest last so README/AUDIT/plot metadata and the complete
    # evidence table are covered by the recorded output hashes.  The plotter
    # refreshes this same manifest after adding figures.
    files = Dict{String, String}()
    for (root, _dirs, names) in walkdir(output_dir)
        for name in names
            name == "manifest.json" && continue
            path = joinpath(root, name)
            files[replace(relpath(path, output_dir), '\\' => '/')] = _sha256_file(path)
        end
    end
    manifest = Dict{String, Any}(
        "schema_version" => SCHEMA_VERSION,
        "verdict" => verdict,
        "source_run_id" => SOURCE_RUN_ID,
        "source_run_conclusion" => data.source_run_conclusion,
        "source_calculation_sha" => expected_calculation_sha,
        "source_postprocess_sha" => expected_postprocess_sha,
        "producer_head_sha" => producer_sha,
        "solver_called" => false,
        "source_job_count" => data.job_count,
        "input_files" => data.input_files,
        "selected_policy" => policy,
        "cost_frontier" => frontier.frontiers,
        "stage_contract" => Dict(
            "stage_a" => Dict("method" => "production_hybrid", "coarse_step" => STAGE_A_COARSE,
                "fine_step" => STAGE_A_FINE, "levels" => [0, 1]),
            "stage_b" => Dict("method" => "independent_oracle", "coarse_step" => STAGE_B_COARSE,
                "fine_step" => STAGE_B_FINE, "levels" => [0]),
            "stage_c_pool" => Dict("method" => "independent_oracle", "step" => STAGE_C_FINE,
                "level" => 1, "disjoint_from_stage_b" => true),
        ),
        "files" => files,
        "oracle_labels_used_for_routing" => false,
        "historical_v1_immutable" => true,
    )
    open(joinpath(output_dir, "manifest.json"), "w") do io
        JSON3.pretty(io, manifest); write(io, '\n')
    end
    manifest
end

function _git_head()
    try
        String(readchomp(`git -C $PROJECT_ROOT rev-parse HEAD`))
    catch
        "unknown"
    end
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
    values
end

function main(args=ARGS)
    options = _parse_args(args)
    expected_calculation_sha = get(options, "expected-calculation-sha", SOURCE_CALCULATION_SHA)
    expected_postprocess_sha = get(options, "expected-source-postprocess-sha", SOURCE_POSTPROCESS_SHA)
    source_run_conclusion = get(options, "source-run-conclusion", "failure")
    producer_sha = get(options, "producer-head-sha", _git_head())
    output_dir = get(options, "output-dir", joinpath(PROJECT_ROOT, "aggregate"))
    data = _load_source(abspath(options["input-dir"]), expected_calculation_sha,
        expected_postprocess_sha, source_run_conclusion)
    frontier = _frontier(data)
    verdict = _verdict(frontier)
    policy = _selected_policy(frontier, verdict)
    manifest = _write_outputs(abspath(output_dir), data, frontier, verdict, policy,
        expected_calculation_sha, expected_postprocess_sha, String(producer_sha))
    println(JSON3.write(Dict("verdict" => verdict, "selected_policy" => policy,
        "source_job_count" => data.job_count, "solver_called" => false,
        "manifest_sha256" => _sha256_file(joinpath(output_dir, "manifest.json")))))
    verdict == "feasible_candidate" ? 0 : 2
end

end # module StageCDensityCertificateFeasibilityV2

if abspath(PROGRAM_FILE) == @__FILE__
    exit(StageCDensityCertificateFeasibilityV2.main())
end
