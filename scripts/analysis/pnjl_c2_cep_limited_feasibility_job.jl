#!/usr/bin/env julia

"""Run one frozen-ξ CEP limited-feasibility fine-pool shard.

The shard is numerical only: it materializes the complete 0.003125 rho pool
at the midpoint of the frozen C2 bracket.  For ξ=0.225 one additional
midpoint is selected from the first hybrid classification.  The oracle labels
are never used to choose a temperature.  Route and oracle gates are replayed
by the solver-free aggregate evaluator.
"""

using Pkg
const PROJECT_ROOT = normpath(pwd())
Pkg.activate(PROJECT_ROOT)

using CSV
using JSON3
using SHA

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
using Main.Models

const SCHEMA_VERSION = "pnjl_c2_cep_limited_feasibility_job_v2"
const CALCULATION_SHA_RE = r"^[0-9a-fA-F]{40}$"
const RHO_FINE_STEP = 0.003125
const HYBRID_LOCAL_STEP = RHO_FINE_STEP / 2
const RHO_MAX = 4.0
# CSV production-eval files use fixed decimal formatting for coordinates.  This
# tolerance covers representation round-off only; it is not a physical gate.
const CSV_COORD_ATOL = 1e-6
const TARGET_XI = (-0.45, -0.34375, -0.275, -0.21875, 0.025, 0.125,
    0.15, 0.2, 0.225, 0.35, 0.3625, 0.38125, 0.39375, 0.4, 0.4375, 0.45, 0.5)
const BRACKET_FILE = joinpath(PROJECT_ROOT, "docs", "analysis", "pnjl",
    "c2_limited_feasibility_v1", "cep_failures.csv")

@inline function _arg(args, name, default=nothing)
    index = findfirst(==(name), args)
    index === nothing && return default
    index == length(args) && throw(ArgumentError("missing value for $(name)"))
    args[index + 1]
end

@inline function _safe_float(value, default=NaN)
    value === nothing || value === missing ? default : try
        result = Float64(value)
        isfinite(result) ? result : default
    catch
        default
    end
end

@inline function _bool(value)
    value === nothing || value === missing ? false :
        value isa Bool ? value : lowercase(strip(String(value))) in ("true", "1", "yes")
end

@inline function _field(record, name::Symbol, default=nothing)
    record === nothing && return default
    try
        if record isa AbstractDict
            haskey(record, name) && return record[name]
            haskey(record, String(name)) && return record[String(name)]
            return default
        end
        hasproperty(record, name) ? getproperty(record, name) : default
    catch
        default
    end
end

function _brackets()
    isfile(BRACKET_FILE) || throw(ArgumentError("missing frozen CEP bracket file: $(BRACKET_FILE)"))
    rows = collect(CSV.File(BRACKET_FILE))
    length(rows) == 17 || throw(ArgumentError("expected 17 frozen CEP brackets"))
    result = NamedTuple[]
    for row in rows
        xi = _safe_float(_field(row, :xi))
        low = _safe_float(_field(row, :T_bracket_low_MeV))
        high = _safe_float(_field(row, :T_bracket_high_MeV))
        isfinite(xi) && isfinite(low) && isfinite(high) && high > low ||
            throw(ArgumentError("invalid frozen CEP bracket: $(row)"))
        push!(result, (xi=xi, low=low, high=high, midpoint=0.5 * (low + high)))
    end
    sort!(result; by=row -> row.xi)
    Set(row.xi for row in result) == Set(TARGET_XI) ||
        throw(ArgumentError("frozen CEP xi matrix does not match v1 contract"))
    result
end

function _config(args)
    scope = String(_arg(args, "--scope", "cep"))
    scope == "cep" || throw(ArgumentError("CEP fine-pool jobs require --scope cep"))
    xi = parse(Float64, String(_arg(args, "--xi", "0.0")))
    xi in TARGET_XI || throw(ArgumentError("unsupported CEP xi $(xi)"))
    calculation_sha = String(_arg(args, "--calculation-sha", ""))
    occursin(CALCULATION_SHA_RE, calculation_sha) ||
        throw(ArgumentError("calculation-sha must be a 40-character SHA"))
    (; scope, xi, calculation_sha,
        output_dir=abspath(String(_arg(args, "--output-dir", joinpath(pwd(), "c2_cep_limited_job")))),
        tag=String(_arg(args, "--tag", "issue130_c2_cep_limited_feasibility_20260812")),
        postprocess_sha=String(_arg(args, "--postprocess-sha", "unknown")),
        source_run_id=String(_arg(args, "--source-run-id", "pending")))
end

@inline function _rho_grid()
    count = Int(round(RHO_MAX / RHO_FINE_STEP))
    Float64.(collect(range(0.0, RHO_MAX; length=count + 1)))
end

function _materialization_row(row, xi, T, calculation_sha, source_path)
    xi_value = _safe_float(_field(row, :xi))
    T_value = _safe_float(_field(row, :T_MeV))
    rho = _safe_float(_field(row, :rho))
    mu = _safe_float(_field(row, :mu_avg_MeV, _field(row, :muq_MeV)))
    residual = _safe_float(_field(row, :residual_norm))
    iterations = _safe_float(_field(row, :iterations))
    converged = _bool(_field(row, :converged, false))
    all(isfinite, (xi_value, T_value, rho, mu, residual, iterations)) && converged ||
        error("invalid production-eval row in $(source_path)")
    isapprox(xi_value, xi; atol=CSV_COORD_ATOL, rtol=0.0) ||
        error("production-eval xi mismatch in $(source_path): $(xi_value) != $(xi)")
    isapprox(T_value, T; atol=CSV_COORD_ATOL, rtol=0.0) ||
        error("production-eval temperature mismatch in $(source_path): $(T_value) != $(T)")
    rho_index = round(Int, rho / RHO_FINE_STEP)
    0 <= rho_index <= round(Int, RHO_MAX / RHO_FINE_STEP) &&
        isapprox(rho, rho_index * RHO_FINE_STEP; atol=3e-7, rtol=0.0) ||
        error("production-eval rho is off the frozen grid in $(source_path): $(rho)")
    (
        # Normalize accepted CSV coordinates before writing the materialized
        # file, so a second parser cannot reintroduce decimal truncation.
        xi=xi, T_MeV=T, rho=rho, muq_MeV=mu,
        residual_norm=residual, iterations=Int(round(iterations)),
        converged=true, finite=true,
        sampling_role="production_eval_materialization_v1", rho_level=0,
        calculation_sha=calculation_sha, rho_index=rho_index,
    )
end

function _materialization_source_files(oracle_dir)
    eval_dir = joinpath(oracle_dir, "production_eval")
    isdir(eval_dir) || error("missing production_eval directory $(eval_dir)")
    paths = filter(path -> isfile(path) && endswith(lowercase(path), ".csv") &&
        occursin("memoized", lowercase(basename(path))), readdir(eval_dir; join=true))
    sort!(paths)
    isempty(paths) && error("missing memoized production-eval files in $(eval_dir)")
    paths
end

function _aggregate_row_count(path)
    isfile(path) || return (rows=0, parse_error="missing aggregate")
    try
        return (rows=length(collect(CSV.File(path))), parse_error=nothing)
    catch error_value
        # The aggregate is not used for reconstruction; retain its state for audit.
        return (rows=max(countlines(path) - 1, 0),
            parse_error=sprint(showerror, error_value))
    end
end

function _materialize_fine_pool(oracle_dir, xi, T, calculation_sha)
    aggregate_path = joinpath(oracle_dir, "trho_scan.csv")
    aggregate = _aggregate_row_count(aggregate_path)
    source_paths = _materialization_source_files(oracle_dir)
    by_key = Dict{Tuple{Float64, Float64, Int}, NamedTuple}()
    source_metadata = Any[]
    for source_path in source_paths
        source_rows = 0
        for row in CSV.File(source_path)
            source_rows += 1
            item = _materialization_row(row, xi, T, calculation_sha, source_path)
            key = (round(item.xi; digits=8), round(item.T_MeV; digits=8), item.rho_index)
            haskey(by_key, key) && error("duplicate production-eval key $(key)")
            by_key[key] = item
        end
        push!(source_metadata, Dict(
            "path" => replace(relpath(source_path, oracle_dir), '\\' => '/'),
            "sha256" => _sha(source_path), "rows" => source_rows,
        ))
    end
    expected_count = length(_rho_grid())
    length(by_key) == expected_count || error(
        "incomplete materialized fine pool at ξ=$(xi), T=$(T): " *
        "$(length(by_key)) != $(expected_count)",
    )
    expected_keys = Set((round(xi; digits=8), round(T; digits=8), index)
        for index in 0:round(Int, RHO_MAX / RHO_FINE_STEP))
    Set(keys(by_key)) == expected_keys || error("materialized fine-pool keys do not match frozen grid")
    rows = collect(values(by_key))
    sort!(rows; by=row -> row.rho)
    output_rows = [NamedTuple{(:xi, :T_MeV, :rho, :muq_MeV, :residual_norm,
        :iterations, :converged, :finite, :sampling_role, :rho_level, :calculation_sha)}(
            (row.xi, row.T_MeV, row.rho, row.muq_MeV, row.residual_norm,
             row.iterations, row.converged, row.finite, row.sampling_role,
             row.rho_level, row.calculation_sha)) for row in rows]
    materialized_path = joinpath(oracle_dir, "trho_scan_materialized.csv")
    CSV.write(materialized_path, output_rows)
    provenance_path = joinpath(oracle_dir, "trho_scan_materialization.json")
    provenance = Dict(
        "schema_version" => "pnjl_c2_cep_fine_pool_materialization_v1",
        "method" => "production_eval_materialization_v1", "solver_called" => false,
        "calculation_sha" => calculation_sha, "xi" => xi, "T_MeV" => T,
        "rho_fine_step" => RHO_FINE_STEP, "rho_max" => RHO_MAX,
        "expected_rows" => expected_count, "materialized_rows" => length(output_rows),
        "recovered_rows" => length(output_rows), "source_aggregate_path" => "trho_scan.csv",
        "source_aggregate_rows" => aggregate.rows,
        "source_aggregate_parse_error" => aggregate.parse_error,
        "source_production_eval_files" => source_metadata,
        "materialized_sha256" => _sha(materialized_path),
    )
    open(provenance_path, "w") do io
        JSON3.pretty(io, provenance)
        write(io, '\n')
    end
    (
        path=materialized_path, provenance_path=provenance_path,
        sha256=_sha(materialized_path), rows=length(output_rows),
        recovered_rows=length(output_rows), aggregate_rows=aggregate.rows,
        aggregate_parse_error=aggregate.parse_error, source_files=source_metadata,
    )
end

function _curve_rows(path, xi, T, calculation_sha)
    isfile(path) || error("missing trho scan $(path)")
    by_rho = Dict{Float64, NamedTuple}()
    for row in CSV.File(path)
        rho = _safe_float(_field(row, :rho))
        mu = _safe_float(_field(row, :mu_avg_MeV, _field(row, :muq_MeV)))
        T_row = _safe_float(_field(row, :T_MeV, T))
        xi_row = _safe_float(_field(row, :xi, xi))
        residual = _safe_float(_field(row, :residual_norm), NaN)
        converged = _bool(_field(row, :converged, false))
        finite = converged && all(isfinite, (xi_row, T_row, rho, mu, residual))
        finite || error("non-finite or non-converged point at ξ=$(xi), T=$(T), rho=$(rho)")
        isapprox(xi_row, xi; atol=CSV_COORD_ATOL, rtol=0.0) ||
            error("trho scan xi mismatch at ξ=$(xi), T=$(T): $(xi_row)")
        isapprox(T_row, T; atol=CSV_COORD_ATOL, rtol=0.0) ||
            error("trho scan temperature mismatch at ξ=$(xi), T=$(T): $(T_row)")
        rho_index = round(Int, rho / RHO_FINE_STEP)
        0 <= rho_index <= round(Int, RHO_MAX / RHO_FINE_STEP) &&
            isapprox(rho, rho_index * RHO_FINE_STEP; atol=3e-7, rtol=0.0) ||
            error("trho scan rho is off the frozen grid at ξ=$(xi), T=$(T): $(rho)")
        item = (
            xi=xi_row, T_MeV=T_row, rho=rho, muq_MeV=mu,
            residual_norm=residual,
            iterations=Int(round(_safe_float(_field(row, :iterations, 0), 0))),
            converged=converged, finite=finite,
            sampling_role="uniform_nested_fine_pool", rho_level=0,
            calculation_sha=calculation_sha,
        )
        haskey(by_rho, rho) && error("duplicate rho key at ξ=$(xi), T=$(T), rho=$(rho)")
        by_rho[rho] = item
    end
    rows = collect(values(by_rho))
    length(rows) == length(_rho_grid()) || error("incomplete fine pool at ξ=$(xi), T=$(T)")
    expected = Set(round(Int, rho / RHO_FINE_STEP) for rho in _rho_grid())
    actual = Set(round(Int, row.rho / RHO_FINE_STEP) for row in rows)
    actual == expected || error("fine pool keys do not match frozen grid at ξ=$(xi), T=$(T)")
    sort!(rows; by=row -> row.rho)
    rows
end

function _stage_b_status(rows)
    fine = [row for row in rows if isapprox(row.rho / 0.00625,
        round(row.rho / 0.00625); atol=3e-7, rtol=0.0)]
    sort!(fine; by=row -> row.rho)
    length(fine) >= 6 || return "ambiguous_near_critical"
    result = Models._classify_s_curve([row.muq_MeV for row in fine], [row.rho for row in fine];
        area_tol_good=1e-4, area_tol_bad=5e-4)
    result.status == :valid && return "confirmed_first_order"
    result.status == :invalid && result.reason == "no_s_shape" && return "confirmed_monotone"
    "ambiguous_near_critical"
end

function _hybrid_status(result)
    records = get(result.diagnostics, "sweep_records", Any[])
    isempty(records) && return (status="ambiguous_near_critical", reason="missing_sweep_record")
    record = first(records)
    status = String(_field(record, :status, "ambiguous_near_critical"))
    reason = String(_field(record, :reason, "unknown"))
    (status=status, reason=reason)
end

function _cache_stats(result)
    diagnostics = get(result.diagnostics, "rho_support_cache", nothing)
    diagnostics === nothing && (diagnostics = get(result.diagnostics, :rho_support_cache, nothing))
    diagnostics === nothing && error("missing rho_support_cache telemetry")
    point_requests = Int(_field(diagnostics, :point_requests, 0))
    cache_hits = Int(_field(diagnostics, :cache_hits, 0))
    unique_solves = Int(_field(diagnostics, :unique_solves, 0))
    failed_points = Int(_field(diagnostics, :failed_points, 0))
    point_requests == unique_solves + cache_hits || error("rho-support cache cost does not reconcile")
    (point_requests=point_requests, cache_hits=cache_hits,
        unique_solves=unique_solves, failed_points=failed_points)
end

function _run_slice(cfg, T, slice_dir)
    mkpath(slice_dir)
    hybrid_telemetry = Models.SolverWorkTelemetry()
    hybrid_dir = joinpath(slice_dir, "hybrid")
    hybrid_started = time_ns()
    hybrid_result = Models.run_production_phase_pipeline(
        :PNJL;
        T_start=T, T_end=T, dT=1.0, rho_grid=Float64.(collect(0.0:0.05:RHO_MAX)), xi=cfg.xi,
        output_dir=hybrid_dir,
        run_id="$(cfg.tag)_hybrid_xi$(cfg.xi)_T$(T)",
        profile=:cep_hybrid_production_shadow, solver_backend=:models,
        reverse_rho=true, seed_policy=:candidates, p_num=24, t_num=8,
        thermo_quadrature_policy=:rs_reduced_adaptive,
        thermo_quadrature_rtol=1e-8, thermo_quadrature_atol=1e-10,
        thermo_quadrature_maxevals=10^7, iterations=80,
        compute_crossover=false, cep_tol=0.125,
        temperature_resolution_target_MeV=0.125, cep_max_bisect_iter=20,
        cep_max_refine_level_rho=4, rho_geometry_convergence=true,
        rho_position_tol_MeV=0.025, rho_density_tol=0.0025,
        rho_maxwell_area_tol=5e-5, rho_refinement_policy=:rho_support_hybrid,
        rho_support_fine_step=RHO_FINE_STEP, rho_support_targeted_cap=12,
        rho_support_config=Models.RhoSupportConfig(),
        rho_hybrid_verification=Models.RhoHybridVerificationConfig(
            local_step=HYBRID_LOCAL_STEP,
            point_ranking_version=:stage_b_features_v1,
            endpoint_policy=:three_crossing_endpoint_local_v2),
        work_telemetry=hybrid_telemetry, memoize_uniform=false, promote_reference=false,
    )
    hybrid_elapsed = (time_ns() - hybrid_started) / 1e9
    hybrid_cache = _cache_stats(hybrid_result)
    oracle_telemetry = Models.SolverWorkTelemetry()
    oracle_dir = joinpath(slice_dir, "oracle")
    oracle_started = time_ns()
    oracle_result = Models.run_production_phase_pipeline(
        :PNJL;
        T_start=T, T_end=T, dT=1.0, rho_grid=_rho_grid(), xi=cfg.xi,
        output_dir=oracle_dir,
        run_id="$(cfg.tag)_oracle_xi$(cfg.xi)_T$(T)",
        profile=:cep_hybrid_production_shadow, solver_backend=:models,
        reverse_rho=true, seed_policy=:candidates, p_num=24, t_num=8,
        thermo_quadrature_policy=:rs_reduced_adaptive,
        thermo_quadrature_rtol=1e-8, thermo_quadrature_atol=1e-10,
        thermo_quadrature_maxevals=10^7, iterations=80,
        compute_crossover=false, cep_tol=0.125,
        temperature_resolution_target_MeV=0.125, cep_max_bisect_iter=20,
        cep_max_refine_level_rho=0, rho_geometry_convergence=false,
        rho_position_tol_MeV=0.025, rho_density_tol=0.0025,
        rho_maxwell_area_tol=5e-5, rho_refinement_policy=:uniform_nested,
        rho_support_fine_step=RHO_FINE_STEP, rho_support_targeted_cap=12,
        rho_support_config=Models.RhoSupportConfig(),
        rho_hybrid_verification=Models.RhoHybridVerificationConfig(),
        work_telemetry=oracle_telemetry, memoize_uniform=true, promote_reference=false,
    )
    oracle_elapsed = (time_ns() - oracle_started) / 1e9
    materialization = _materialize_fine_pool(oracle_dir, cfg.xi, T, cfg.calculation_sha)
    oracle_rows = _curve_rows(materialization.path, cfg.xi, T, cfg.calculation_sha)
    oracle_cache = _cache_stats(oracle_result)
    costs = [
        (method="hybrid", cache=hybrid_cache, telemetry=Models.solver_work_snapshot(hybrid_telemetry),
            runner_seconds=hybrid_elapsed),
        (method="oracle", cache=oracle_cache, telemetry=Models.solver_work_snapshot(oracle_telemetry),
            runner_seconds=oracle_elapsed),
    ]
    all(item -> item.cache.failed_points == 0, costs) || error("solver failure recorded in CEP slice")
    oracle_rows, _hybrid_status(hybrid_result), costs, materialization
end

function _sha(path)
    bytes2hex(open(sha256, path))
end

function _write_csv(path, rows)
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _run_job(cfg)
    brackets = _brackets()
    bracket = first(row for row in brackets if row.xi == cfg.xi)
    mkpath(cfg.output_dir)
    started = time_ns()
    plan = [(T=bracket.midpoint, role="midpoint", low=bracket.low, high=bracket.high,
        selection_reason="frozen_bracket_midpoint")]
    curve_rows = NamedTuple[]
    slice_rows = NamedTuple[]
    cost_rows = NamedTuple[]
    materialization_rows = Any[]
    index = 1
    while index <= length(plan)
        item = plan[index]
        rows, hybrid, costs, materialization = _run_slice(cfg, item.T,
            joinpath(cfg.output_dir, "slices", "T_$(replace(string(item.T), "." => "p"))"))
        push!(materialization_rows, Dict(
            "xi" => cfg.xi, "T_MeV" => item.T,
            "path" => replace(relpath(materialization.path, cfg.output_dir), '\\' => '/'),
            "provenance_path" => replace(relpath(materialization.provenance_path, cfg.output_dir), '\\' => '/'),
            "sha256" => materialization.sha256, "rows" => materialization.rows,
            "recovered_rows" => materialization.recovered_rows,
            "aggregate_rows" => materialization.aggregate_rows,
            "aggregate_parse_error" => materialization.aggregate_parse_error,
            "source_files" => materialization.source_files,
        ))
        status = hybrid.status
        push!(curve_rows, rows...)
        for item_cost in costs
            snapshot = item_cost.telemetry
            push!(cost_rows, (
                xi=cfg.xi, T_MeV=item.T, method=item_cost.method,
                calculation_sha=cfg.calculation_sha,
                unique_solves=item_cost.cache.unique_solves,
                point_requests=item_cost.cache.point_requests,
                requested_point_calls=item_cost.cache.point_requests,
                cache_hits=item_cost.cache.cache_hits,
                failed_points=item_cost.cache.failed_points,
                equilibrium_requests=snapshot.equilibrium_requests,
                fixedrho_requests=snapshot.fixedrho_requests,
                residual_calls=snapshot.nlsolve_f_calls + snapshot.postprocess_residual_calls,
                jacobian_calls=snapshot.nlsolve_g_calls,
                newton_iterations=snapshot.newton_iterations,
                runner_seconds=item_cost.runner_seconds,
                fallback_count=snapshot.root_fallbacks,
                retry_count=snapshot.scan_retries,
            ))
        end
        push!(slice_rows, (xi=cfg.xi, T_MeV=item.T, target_role=item.role,
            bracket_low_MeV=item.low, bracket_high_MeV=item.high,
            hybrid_status=status, hybrid_reason=hybrid.reason,
            selected_by="hybrid_production_status", selection_reason=item.selection_reason,
            oracle_labels_used_for_routing=false,
            oracle_materialization_method="production_eval_materialization_v1",
            oracle_materialized_rows=materialization.rows,
            oracle_recovered_rows=materialization.recovered_rows,
            oracle_aggregate_rows=materialization.aggregate_rows,
            oracle_materialized_sha256=materialization.sha256,
            oracle_materialization_provenance=replace(
                relpath(materialization.provenance_path, cfg.output_dir), '\\' => '/'),
            rho_rows=length(rows),
            finite_and_converged=true, slice_index=index))
        if cfg.xi == 0.225 && index == 1 &&
                status in ("confirmed_first_order", "confirmed_monotone")
            extra = status == "confirmed_first_order" ?
                0.5 * (item.T + bracket.high) : 0.5 * (bracket.low + item.T)
            extra_low = status == "confirmed_first_order" ? item.T : bracket.low
            extra_high = status == "confirmed_first_order" ? bracket.high : item.T
            push!(plan, (T=extra, role="extra_midpoint", low=extra_low, high=extra_high,
                selection_reason="hybrid_midpoint_status_$(status)"))
        end
        index += 1
    end
    sort!(curve_rows; by=row -> (row.T_MeV, row.rho))
    elapsed = (time_ns() - started) / 1e9
    _write_csv(joinpath(cfg.output_dir, "fine_pool.csv"), curve_rows)
    _write_csv(joinpath(cfg.output_dir, "slice_metrics.csv"), slice_rows)
    all(item.failed_points == 0 for item in cost_rows) || error("CEP shard has failed solver points")
    all(item.point_requests == item.unique_solves + item.cache_hits for item in cost_rows) ||
        error("CEP shard cache cost does not reconcile")
    _write_csv(joinpath(cfg.output_dir, "method_costs.csv"), cost_rows)
    files = Dict(name => _sha(joinpath(cfg.output_dir, name)) for name in
        ("fine_pool.csv", "slice_metrics.csv", "method_costs.csv"))
    summary = Dict(
        "schema_version" => SCHEMA_VERSION, "scope" => cfg.scope, "xi" => cfg.xi,
        "method" => "hybrid_and_oracle_fine_pool", "tag" => cfg.tag,
        "source_run_id" => cfg.source_run_id, "calculation_sha" => cfg.calculation_sha,
        "postprocess_sha" => cfg.postprocess_sha, "workflow_head_sha" => cfg.postprocess_sha,
        "rho_fine_step" => RHO_FINE_STEP, "rho_max" => RHO_MAX,
        "hybrid_local_step" => HYBRID_LOCAL_STEP,
        "frozen_bracket" => Dict("low_MeV" => bracket.low, "high_MeV" => bracket.high,
            "midpoint_MeV" => bracket.midpoint),
        "temperatures" => [row.T_MeV for row in slice_rows],
        "oracle_materialization" => materialization_rows,
        "solver_called" => true, "finite_and_converged_final" => true,
        "cache" => Dict("methods" => Dict(
            item.method => Dict("unique_solves" => item.unique_solves,
                "point_requests" => item.point_requests, "cache_hits" => item.cache_hits,
                "failed_points" => item.failed_points) for item in cost_rows)),
        "telemetry" => Dict(item.method => Dict(
            string(field) => getproperty(item, field) for field in propertynames(item))
            for item in cost_rows),
        "runner_seconds" => elapsed, "files" => files,
    )
    open(joinpath(cfg.output_dir, "job_summary.json"), "w") do io
        JSON3.write(io, summary); write(io, '\n')
    end
    open(joinpath(cfg.output_dir, "manifest.json"), "w") do io
        JSON3.pretty(io, merge(summary, Dict("provenance" => Dict(
            "calculation_sha" => cfg.calculation_sha, "postprocess_sha" => cfg.postprocess_sha,
            "reference_write" => false, "route_before_oracle_gate" => true,
            "oracle_labels_used_for_routing" => false,
            "hybrid_local_step" => HYBRID_LOCAL_STEP,
            "hybrid_local_step_contract" => "rho_support_fine_step/2"))))
        write(io, '\n')
    end
    println(JSON3.write(summary))
end

if abspath(PROGRAM_FILE) == @__FILE__
    _run_job(_config(ARGS))
end
