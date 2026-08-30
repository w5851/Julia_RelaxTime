#!/usr/bin/env julia

"""Solver-free CEP limited-feasibility replay.

The numerical shard owns equilibrium calls and only materializes frozen
bracket midpoints.  This evaluator reconstructs full fine curves, runs the
calculation SHA's production-parity classifier and Maxwell implementation,
and uses oracle rows only after route selection for the final gate.
"""

using Pkg
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(PROJECT_ROOT)
using CSV
using JSON3
using SHA
using Statistics

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
const MODELS = Main.Models

const SCHEMA_VERSION = "pnjl_c2_cep_limited_feasibility_v2"
const JOB_SCHEMA = "pnjl_c2_cep_limited_feasibility_job_v2"
const CALCULATION_SHA_RE = r"^[0-9a-f]{40}$"
const RHO_FINE_STEP = 0.003125
const STAGE_B_STEP = 0.00625
const TARGET_XI = (-0.45, -0.34375, -0.275, -0.21875, 0.025, 0.125,
    0.15, 0.2, 0.225, 0.35, 0.3625, 0.38125, 0.39375, 0.4, 0.4375, 0.45, 0.5)
const WIDTH_GATE = 0.1
const COST_STOP_MINUTES = 200.0
const PRODUCTION_CONFIG = MODELS.ProductionPipelineConfig(
    area_tol_good=1e-4,
    area_tol_bad=5e-4,
    rho_geometry_convergence=true,
    rho_position_tol_MeV=0.025,
    rho_density_tol=0.0025,
    rho_maxwell_area_tol=5e-5,
)
const PRODUCTION_MAXWELL_OPTIONS = MODELS._production_maxwell_options(PRODUCTION_CONFIG)

@inline function _arg(args, name, default=nothing)
    index = findfirst(==(name), args)
    index === nothing && return default
    index == length(args) && error("missing value for $(name)")
    args[index + 1]
end

function _parse_target_xi(value::AbstractString)
    isempty(strip(value)) && return TARGET_XI
    values = try
        [parse(Float64, strip(item)) for item in split(value, ',') if !isempty(strip(item))]
    catch
        error("invalid --target-xi value: $(value)")
    end
    isempty(values) && error("--target-xi must contain at least one xi")
    length(unique(values)) == length(values) || error("--target-xi must contain unique xi values")
    all(xi -> xi in TARGET_XI, values) ||
        error("--target-xi contains an xi outside the frozen matrix")
    sort!(values)
    Tuple(values)
end

@inline function _field(row, name::Symbol, default=nothing)
    row === nothing && return default
    try
        if row isa AbstractDict
            haskey(row, name) && return row[name]
            haskey(row, String(name)) && return row[String(name)]
            return default
        end
        hasproperty(row, name) ? getproperty(row, name) : default
    catch
        default
    end
end

@inline _field(row, name::AbstractString, default=nothing) =
    _field(row, Symbol(name), default)

@inline function _float(value, default=NaN)
    value === nothing || value === missing ? default : try
        result = Float64(value)
        isfinite(result) ? result : default
    catch
        default
    end
end

# Keep temperature lookup keys stable across CSV round-off while remaining
# compatible with Julia 1.12, where `round(x, digits)` is not a valid call.
@inline _temperature_key(value) = round(Float64(value); digits=8)

@inline function _bool(value)
    value === nothing || value === missing ? false : value isa Bool ? value :
        lowercase(strip(String(value))) in ("true", "1", "yes")
end

function _sha(path::String)
    bytes2hex(open(sha256, path))
end

function _manifest_paths(input_dir::String)
    paths = String[]
    for (root, _dirs, files) in walkdir(input_dir)
        "manifest.json" in files && push!(paths, joinpath(root, "manifest.json"))
    end
    sort(paths)
end

function _required_file(manifest, manifest_path::String, name::String)
    files = _field(manifest, :files, nothing)
    declared = _field(files, Symbol(name), nothing)
    declared === nothing && error("manifest missing files.$(name): $(manifest_path)")
    path = joinpath(dirname(manifest_path), name)
    isfile(path) || error("missing $(name): $(manifest_path)")
    String(declared) == _sha(path) || error("hash mismatch for $(name): $(manifest_path)")
    path
end

function _read_jobs(input_dir::String, expected_sha::String, expected_postprocess::String,
        target_xi=TARGET_XI)
    paths = _manifest_paths(input_dir)
    length(paths) == length(target_xi) ||
        error("expected $(length(target_xi)) CEP shards, got $(length(paths))")
    jobs = Any[]
    seen = Set{Float64}()
    for path in paths
        manifest = JSON3.read(read(path, String))
        String(_field(manifest, :schema_version, "")) == JOB_SCHEMA || error("unexpected job schema $(path)")
        String(_field(manifest, :scope, "")) == "cep" || error("unexpected job scope $(path)")
        sha = lowercase(String(_field(manifest, :calculation_sha, "")))
        sha == expected_sha || error("calculation SHA mismatch $(path)")
        post = String(_field(manifest, :postprocess_sha, ""))
        isempty(expected_postprocess) || post == expected_postprocess || error("postprocess SHA mismatch $(path)")
        xi = _float(_field(manifest, :xi))
        xi in target_xi && !(xi in seen) || error("duplicate/unexpected xi $(xi)")
        push!(seen, xi)
        pool_path = _required_file(manifest, path, "fine_pool.csv")
        slices_path = _required_file(manifest, path, "slice_metrics.csv")
        costs_path = _required_file(manifest, path, "method_costs.csv")
        _bool(_field(manifest, :solver_called, false)) || error("CEP source must record solver_called=true")
        provenance = _field(manifest, :provenance, nothing)
        _bool(_field(provenance, :route_before_oracle_gate, false)) ||
            error("CEP source must record route_before_oracle_gate=true")
        points = NamedTuple[]
        keys = Set{Tuple{Float64, Float64}}()
        for row in CSV.File(pool_path)
            T = _float(_field(row, :T_MeV)); rho = _float(_field(row, :rho)); mu = _float(_field(row, :muq_MeV))
            finite = _bool(_field(row, :finite, false)) && _bool(_field(row, :converged, false))
            finite && all(isfinite, (T, rho, mu, _float(_field(row, :residual_norm)))) || error("invalid source point $(path)")
            key = (T, rho)
            key in keys && error("duplicate source point $(key)")
            push!(keys, key)
            push!(points, (xi=xi, T_MeV=T, rho=rho, muq_MeV=mu,
                residual_norm=_float(_field(row, :residual_norm)), finite=true, converged=true))
        end
        isempty(points) && error("empty fine pool $(path)")
        slices = Any[collect(CSV.File(slices_path))...]
        isempty(slices) && error("missing CEP slice metrics $(path)")
        length(slices) <= 2 || error("too many CEP slices $(path)")
        all(_bool(_field(row, :oracle_labels_used_for_routing, false)) == false for row in slices) ||
            error("oracle labels must not route CEP slices $(path)")
        costs = Any[collect(CSV.File(costs_path))...]
        push!(jobs, (; xi, points, slices, costs, manifest,
            files=[(path=replace(relpath(p, input_dir), '\\' => '/'), sha256=_sha(p))
                for p in (pool_path, slices_path, costs_path, path)]))
    end
    Set(job.xi for job in jobs) == Set(target_xi) || error("incomplete CEP xi matrix")
    jobs
end

function _curve(points)
    ordered = sort(points; by=row -> row.rho)
    length(ordered) >= 6 || return nothing
    length(unique(row.rho for row in ordered)) == length(ordered) || return nothing
    all(isfinite(row.rho) && isfinite(row.muq_MeV) for row in ordered) || return nothing
    (rho=[row.rho for row in ordered], mu=[row.muq_MeV for row in ordered], points=ordered)
end

function _evaluate(points)
    curve = _curve(points)
    curve === nothing && return (status=:ambiguous_near_critical, reason="curve_too_short",
        result=MODELS.MaxwellResult(), classification=nothing, candidate_count=0, crossing_count=0)
    cres = MODELS._classify_s_curve(curve.mu, curve.rho;
        maxwell_options=PRODUCTION_MAXWELL_OPTIONS, area_tol_good=1e-4, area_tol_bad=5e-4)
    raw_maxwell = get(cres, :maxwell, MODELS.MaxwellResult())
    maxwell = MODELS._production_maxwell_or_empty(cres, PRODUCTION_CONFIG)
    details = raw_maxwell.details
    candidate_count = Int(_field(details, :candidate_count, _field(details, "candidate_count", 0)))
    crossing_count = Int(_field(details, :crossing_count, _field(details, "crossing_count", 0)))
    if cres.status == :invalid && cres.reason == "no_s_shape"
        return (status=:confirmed_monotone, reason="no_s_shape", result=maxwell,
            classification=cres, candidate_count=0, crossing_count=0)
    elseif cres.status == :valid && maxwell.converged && maxwell.mu_transition !== nothing &&
            candidate_count == 1 && crossing_count == 3
        return (status=:confirmed_first_order, reason="valid_three_crossing_candidate",
            result=maxwell, classification=cres, candidate_count=candidate_count,
            crossing_count=crossing_count)
    end
    reason = String(_field(details, :reason, _field(details, "reason", "maxwell_unresolved")))
    (status=:ambiguous_near_critical, reason=reason, result=maxwell,
        classification=cres, candidate_count=candidate_count, crossing_count=crossing_count)
end

function _temperature_rows(job)
    temperatures = sort(unique(row.T_MeV for row in job.points))
    rows = NamedTuple[]
    for T in temperatures
        points = [row for row in job.points if isapprox(row.T_MeV, T; atol=1e-8, rtol=0.0)]
        push!(rows, (xi=job.xi, T_MeV=T, points=points, eval=_evaluate(points)))
    end
    rows
end

@inline function _status(value)
    text = lowercase(strip(String(value)))
    occursin("first_order", text) || text == "valid" ? "confirmed_first_order" :
        (occursin("monotone", text) || text == "no_s_shape" ? "confirmed_monotone" :
        "ambiguous_near_critical")
end

function _frozen_bracket(job)
    manifest_bracket = _field(job.manifest, :frozen_bracket, nothing)
    low = _float(_field(manifest_bracket, :low_MeV, _field(manifest_bracket, :low, NaN)))
    high = _float(_field(manifest_bracket, :high_MeV, _field(manifest_bracket, :high, NaN)))
    if !(isfinite(low) && isfinite(high) && high > low)
        first_slice = first(job.slices)
        low = _float(_field(first_slice, :bracket_low_MeV))
        high = _float(_field(first_slice, :bracket_high_MeV))
    end
    isfinite(low) && isfinite(high) && high > low || error("invalid frozen bracket for xi=$(job.xi)")
    (low=low, high=high)
end

function _endpoint_replay(job, rows)
    bracket = _frozen_bracket(job)
    ordered = sort(job.slices; by=row -> _float(_field(row, :T_MeV)))
    oracle = Dict(_temperature_key(row.T_MeV) => String(row.eval.status) for row in rows)
    hybrid = Dict(_temperature_key(_float(_field(row, :T_MeV))) => _status(_field(row, :hybrid_status,
        "ambiguous_near_critical")) for row in ordered)
    low = bracket.low
    high = bracket.high
    low_status = "confirmed_first_order"
    high_status = "confirmed_monotone"
    low_source = "frozen_c2_bracket_low"
    high_source = "frozen_c2_bracket_high"
    invalid_selection = false
    invalid_reason = ""
    previous_low, previous_high = low, high
    for row in ordered
        T = _float(_field(row, :T_MeV))
        expected_midpoint = 0.5 * (previous_low + previous_high)
        if !isfinite(T) || abs(T - expected_midpoint) > 1e-7
            invalid_selection = true
            invalid_reason = "midpoint_not_inside_current_frozen_bracket"
            break
        end
        # Temperature routing is a property of the hybrid run.  Oracle labels
        # are compared only after the route has been reconstructed.
        status = get(hybrid, _temperature_key(T), "ambiguous_near_critical")
        if status == "confirmed_first_order"
            T > low + 1e-8 || (invalid_selection = true; invalid_reason = "first_order_midpoint_not_above_low")
            low = T; low_status = status; low_source = "oracle_midpoint"
        elseif status == "confirmed_monotone"
            T < high - 1e-8 || (invalid_selection = true; invalid_reason = "monotone_midpoint_not_below_high")
            high = T; high_status = status; high_source = "oracle_midpoint"
        else
            invalid_selection = true
            invalid_reason = "oracle_midpoint_ambiguous"
            break
        end
        previous_low, previous_high = low, high
    end
    hybrid_mismatch = any(get(hybrid, key, "missing") != get(oracle, key, "missing")
        for key in keys(oracle))
    width = high - low
    endpoint_status = invalid_selection ? "ambiguous" :
        (low_status == "confirmed_first_order" && high_status == "confirmed_monotone" &&
         width <= WIDTH_GATE ? "pass" : "unresolved")
    (; bracket_low=bracket.low, bracket_high=bracket.high, low, high,
        width, low_status, high_status, low_source, high_source,
        hybrid_mismatch, invalid_selection, invalid_reason, endpoint_status,
        hybrid_states=hybrid, oracle_states=oracle)
end

function _pipeline_cost(job)
    total = 0.0
    complete = true
    rows = NamedTuple[]
    for row in job.costs
        value = _float(_field(row, :runner_seconds, NaN))
        isfinite(value) || (complete = false)
        total += isfinite(value) ? value : 0.0
        method = String(_field(row, :method, ""))
        method in ("hybrid", "oracle") || (complete = false)
        unique_solves = _float(_field(row, :unique_solves, NaN))
        point_requests = _float(_field(row, :point_requests,
            _field(row, :requested_point_calls, NaN)))
        cache_hits = _float(_field(row, :cache_hits, NaN))
        all(isfinite, (unique_solves, point_requests, cache_hits)) || (complete = false)
        point_requests == unique_solves + cache_hits || (complete = false)
        push!(rows, (method=method, unique_solves=unique_solves,
            point_requests=point_requests, cache_hits=cache_hits,
            runner_seconds=value))
    end
    methods = Set(row.method for row in rows)
    methods == Set(("hybrid", "oracle")) || (complete = false)
    (runner_seconds=total, telemetry_complete=complete, rows=rows,
        unique_solves=sum(row.unique_solves for row in rows),
        point_requests=sum(row.point_requests for row in rows),
        cache_hits=sum(row.cache_hits for row in rows))
end

function _job_result(job)
    rows = _temperature_rows(job)
    endpoints = _endpoint_replay(job, rows)
    costs = _pipeline_cost(job)
    state = endpoints.endpoint_status
    (xi=job.xi, rows=rows, hybrid_states=endpoints.hybrid_states,
        T_last_first_order_MeV=endpoints.low,
        T_first_monotone_MeV=endpoints.high,
        bracket_width_T_MeV=endpoints.width, endpoint_status=state,
        bracket_low_MeV=endpoints.bracket_low, bracket_high_MeV=endpoints.bracket_high,
        low_status=endpoints.low_status, high_status=endpoints.high_status,
        low_source=endpoints.low_source, high_source=endpoints.high_source,
        endpoint_reason=endpoints.invalid_reason,
        hybrid_mismatch=endpoints.hybrid_mismatch,
        candidate_failures=sum(row.eval.candidate_count > 1 ||
            (row.eval.candidate_count > 0 && row.eval.crossing_count != 3) for row in rows),
        finite_and_converged=all(all(row -> row.finite && row.converged, item.points) for item in rows),
        telemetry_complete=costs.telemetry_complete, runner_seconds=costs.runner_seconds,
        unique_solves=costs.unique_solves, point_requests=costs.point_requests,
        cache_hits=costs.cache_hits)
end

function _verdict(results)
    any(!result.finite_and_converged for result in results) && return "cep_solver_or_curve_failure"
    sum(result.runner_seconds for result in results) / 60 > COST_STOP_MINUTES && return "cep_runner_cost_exceeded"
    any(!result.telemetry_complete for result in results) && return "cep_telemetry_incomplete"
    any(result.hybrid_mismatch for result in results) && return "cep_hybrid_integration_failed"
    any(result.candidate_failures > 0 for result in results) && return "cep_maxwell_candidate_inconclusive"
    any(result.endpoint_status == "ambiguous" for result in results) && return "cep_oracle_inconclusive"
    all(result.endpoint_status == "pass" for result in results) || return "cep_ambiguity_inconclusive"
    "cep_feasible_candidate"
end

function _write_csv(path, rows)
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _write_outputs(output_dir::String, jobs, results, verdict, expected_sha,
        expected_postprocess, source_run_id, target_xi, workflow_schema)
    mkpath(output_dir)
    _write_csv(joinpath(output_dir, "cep_bracket_results.csv"), [(
        xi=result.xi, T_last_first_order_MeV=result.T_last_first_order_MeV,
        T_first_monotone_MeV=result.T_first_monotone_MeV,
        bracket_width_T_MeV=result.bracket_width_T_MeV,
        endpoint_status=result.endpoint_status,
        hybrid_mismatch=result.hybrid_mismatch,
        finite_and_converged=result.finite_and_converged,
        runner_seconds=result.runner_seconds, unique_solves=result.unique_solves,
        point_requests=result.point_requests, cache_hits=result.cache_hits,
        telemetry_complete=result.telemetry_complete) for result in results])
    _write_csv(joinpath(output_dir, "temperature_states.csv"), [(
        xi=result.xi, T_MeV=row.T_MeV, status=String(row.eval.status), reason=row.eval.reason,
        point_count=length(row.points), hybrid_status=get(result.hybrid_states, _temperature_key(row.T_MeV), "missing"))
        for result in results for row in result.rows])
    _write_csv(joinpath(output_dir, "method_costs.csv"), [(
        xi=result.xi, method="aggregate", unique_solves=result.unique_solves,
        point_requests=result.point_requests, cache_hits=result.cache_hits,
        runner_seconds=result.runner_seconds, telemetry_complete=result.telemetry_complete)
        for result in results])
    _write_csv(joinpath(output_dir, "cost_frontier.csv"), [(
        scope="cep", cap=12, bracket_count=length(results),
        runner_minutes=sum(result.runner_seconds for result in results) / 60,
        cost_stop_minutes=200, cost_gate=sum(result.runner_seconds for result in results) / 60 <= COST_STOP_MINUTES,
        unique_solves=sum(result.unique_solves for result in results),
        point_requests=sum(result.point_requests for result in results),
        cache_hits=sum(result.cache_hits for result in results),
        point_request_reconciliation=sum(result.point_requests for result in results) ==
            sum(result.unique_solves for result in results) + sum(result.cache_hits for result in results),
        verdict=verdict)])
    _write_csv(joinpath(output_dir, "route_trace.csv"), [(
        xi=result.xi, T_MeV=row.T_MeV, route="stage_b_features_v1",
        source="full_fine_pool", oracle_gate="post_route_only",
        status=String(row.eval.status), reason=row.eval.reason) for result in results for row in result.rows])
    write(joinpath(output_dir, "README.md"), "# C2 CEP limited feasibility v2\n\n" *
        "This is diagnostic-only evidence. Frozen bracket midpoints were selected before oracle evaluation; " *
        "the aggregate replay uses the calculation SHA classifier and Maxwell implementation.\n")
    write(joinpath(output_dir, "AUDIT.md"), "# CEP limited feasibility audit\n\n" *
        "No reference, production or transport artifact is modified.\n")
    claims = [
        (claim_id="cep_bracket_gate", status=verdict == "cep_feasible_candidate" ? "pass" : "inconclusive", evidence="cep_bracket_results.csv"),
        (claim_id="cep_classification", status=verdict == "cep_feasible_candidate" ? "pass" : "inconclusive", evidence="temperature_states.csv"),
        (claim_id="production_promotion", status="not_claimed", evidence="author review required"),
    ]
    _write_csv(joinpath(output_dir, "claim_ledger.csv"), claims)
    files = Dict{String, String}()
    for (root, _dirs, names) in walkdir(output_dir)
        for name in names
            name == "manifest.json" && continue
            path = joinpath(root, name)
            files[replace(relpath(path, output_dir), '\\' => '/')] = _sha(path)
        end
    end
    manifest = Dict(
        "schema_version" => SCHEMA_VERSION, "workflow_schema" => workflow_schema,
        "scope" => "cep", "verdict" => verdict,
        "source_run_id" => source_run_id, "source_job_count" => length(jobs),
        "source_calculation_sha" => expected_sha, "source_postprocess_sha" => expected_postprocess,
        "solver_called" => false, "oracle_labels_used_for_routing" => false,
        "frozen_bracket_count" => length(target_xi),
        "max_new_temperature_slices" => 2 * length(target_xi),
        "target_xi" => collect(target_xi),
        "node_step_MeV" => 0.0625, "window_extension_MeV" => 0.25,
        "width_gate_MeV" => WIDTH_GATE, "cost_stop_runner_minutes" => COST_STOP_MINUTES,
        "input_files" => reduce(vcat, [job.files for job in jobs]; init=Any[]), "files" => files)
    open(joinpath(output_dir, "manifest.json"), "w") do io
        JSON3.pretty(io, manifest); write(io, '\n')
    end
    manifest
end

function main(args=ARGS)
    options = Dict{String, String}()
    index = 1
    while index <= length(args)
        arg = args[index]
        startswith(arg, "--") || error("unexpected argument $(arg)")
        key = replace(first(split(arg, "="; limit=2)), "--" => "")
        if occursin("=", arg)
            options[key] = split(arg, "="; limit=2)[2]
        else
            index == length(args) && error("missing value for --$(key)")
            options[key] = args[index + 1]; index += 1
        end
        index += 1
    end
    input_dir = get(options, "input-dir", "")
    isempty(input_dir) && error("--input-dir is required")
    expected_sha = lowercase(get(options, "expected-calculation-sha", ""))
    occursin(CALCULATION_SHA_RE, expected_sha) || error("expected calculation SHA must be lowercase 40-hex")
    source_run_id = get(options, "source-run-id", "")
    occursin(r"^\d+$", source_run_id) || error("source-run-id must be numeric")
    target_xi = _parse_target_xi(get(options, "target-xi", ""))
    workflow_schema = get(options, "workflow-schema", SCHEMA_VERSION)
    jobs = _read_jobs(abspath(input_dir), expected_sha,
        get(options, "expected-source-postprocess-sha", ""), target_xi)
    results = [_job_result(job) for job in jobs]
    verdict = _verdict(results)
    output_dir = abspath(get(options, "output-dir", joinpath(PROJECT_ROOT, "c2_cep_limited_aggregate")))
    manifest = _write_outputs(output_dir, jobs, results, verdict, expected_sha,
        get(options, "expected-source-postprocess-sha", ""), source_run_id, target_xi,
        workflow_schema)
    println(JSON3.write(Dict("verdict" => verdict, "source_job_count" => length(jobs),
        "solver_called" => false, "manifest_sha256" => _sha(joinpath(output_dir, "manifest.json")))))
    verdict == "cep_feasible_candidate" ? 0 : 2
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
