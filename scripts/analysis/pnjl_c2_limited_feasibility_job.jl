#!/usr/bin/env julia

"""Run one ξ shard for the Issue #130 C2 limited-feasibility workflow.

The calculation checkout is immutable and the shard produces one complete
0.003125 rho pool for every requested temperature at that ξ.  Selection and
certificate decisions are deliberately left to the solver-free evaluator.
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

const SCHEMA_VERSION = "pnjl_c2_limited_feasibility_job_v1"
const CALCULATION_SHA_RE = r"^[0-9a-fA-F]{40}$"
const RHO_FINE_STEP = 0.003125
const RHO_MAX = 4.0
const METHOD = "uniform_nested_fine_pool"

const DENSITY_ANCHORS = Dict(
    -0.5 => [60.0, 160.0],
    -0.35 => [51.0],
    -0.25 => [41.0],
    -0.2 => [41.0],
    -0.15 => [41.0],
    -0.1 => [41.0],
    0.0 => [51.0, 60.0, 145.0],
    0.3 => [21.0],
    0.35 => [51.0, 101.0],
    0.5 => [60.0, 120.0],
)

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

@inline function _safe_int(value, default=0)
    result = _safe_float(value, NaN)
    isfinite(result) ? Int(round(result)) : default
end

@inline function _bool(value)
    (value === nothing || value === missing) && return false
    value isa Bool && return value
    lowercase(strip(String(value))) in ("true", "1", "yes")
end

@inline function _field(record, name::Symbol, default=nothing)
    record === nothing && return default
    try
        hasproperty(record, name) ? getproperty(record, name) : default
    catch
        default
    end
end

function _config(args)
    scope = String(_arg(args, "--scope", "density"))
    scope == "density" || throw(ArgumentError("numerical fine-pool jobs currently require --scope density"))
    xi = parse(Float64, String(_arg(args, "--xi", "0.0")))
    haskey(DENSITY_ANCHORS, xi) || throw(ArgumentError("unsupported density ξ $(xi)"))
    calculation_sha = String(_arg(args, "--calculation-sha", ""))
    occursin(CALCULATION_SHA_RE, calculation_sha) ||
        throw(ArgumentError("calculation-sha must be a 40-character SHA"))
    output_dir = abspath(String(_arg(args, "--output-dir", joinpath(pwd(), "c2_limited_job"))))
    tag = String(_arg(args, "--tag", "issue130_c2_density_feasibility_20260812"))
    postprocess_sha = String(_arg(args, "--postprocess-sha", "unknown"))
    source_run_id = String(_arg(args, "--source-run-id", "pending"))
    (; scope, xi, calculation_sha, output_dir, tag, postprocess_sha, source_run_id)
end

function _rho_grid()
    count = Int(round(RHO_MAX / RHO_FINE_STEP))
    Float64.(collect(range(0.0, RHO_MAX; length=count + 1)))
end

function _curve_rows(path, xi, T, calculation_sha)
    isfile(path) || error("missing trho scan $(path)")
    rows = NamedTuple[]
    seen = Set{Tuple{Float64, Float64, Float64}}()
    for row in CSV.File(path)
        rho = _safe_float(_field(row, :rho))
        mu = _safe_float(_field(row, :mu_avg_MeV, _field(row, :muq_MeV)))
        T_row = _safe_float(_field(row, :T_MeV, T))
        xi_row = _safe_float(_field(row, :xi, xi))
        residual = _safe_float(_field(row, :residual_norm, NaN))
        key = (xi_row, T_row, rho)
        key in seen && error("duplicate rho key $(key)")
        push!(seen, key)
        converged = _bool(_field(row, :converged, false))
        finite = converged && isfinite(rho) && isfinite(mu) && isfinite(residual)
        finite || error("non-finite or non-converged point at ξ=$(xi), T=$(T), rho=$(rho)")
        push!(rows, (
            xi=xi_row, T_MeV=T_row, rho=rho, muq_MeV=mu,
            residual_norm=residual, iterations=Int(round(_safe_float(_field(row, :iterations, 0), 0))),
            converged=converged, finite=finite, sampling_role="uniform_nested_fine_pool",
            rho_level=0, calculation_sha=calculation_sha,
        ))
    end
    sort!(rows; by=row -> (row.T_MeV, row.rho))
    rows
end

function _slice_row(result, xi, T, calculation_sha)
    diagnostics = result.diagnostics
    records = get(diagnostics, "sweep_records", NamedTuple[])
    record = findfirst(row -> isapprox(_safe_float(_field(row, :T_MeV)), T; atol=1e-8, rtol=0.0), records)
    record = record === nothing ? nothing : records[record]
    failure = Int(_safe_float(_field(record, :stats_failure, get(diagnostics, "scan_failure", 0)), 0))
    cache = get(diagnostics, "rho_support_cache", Dict{String, Any}())
    cache = cache isa AbstractDict ? cache : Dict{String, Any}()
    (
        xi=xi, method=METHOD, T_MeV=T, calculation_sha=calculation_sha,
        result_status=String(_field(record, :status, result.cep.result_status)),
        stage_status=String(_field(record, :stage_used, :uniform_nested)),
        rho_hadron=_safe_float(_field(record, :rho_hadron)),
        rho_quark=_safe_float(_field(record, :rho_quark)),
        rho_spinodal_hadron=_safe_float(_field(record, :rho_spinodal_hadron)),
        rho_spinodal_quark=_safe_float(_field(record, :rho_spinodal_quark)),
        mu_transition_MeV=_safe_float(_field(record, :mu_transition_MeV, _field(record, :mu_transition))),
        geometry_converged=_bool(_field(record, :geometry_converged, false)),
        position_error_MeV=_safe_float(_field(record, :position_error_MeV)),
        density_error=_safe_float(_field(record, :density_error)),
        maxwell_area_gate=_safe_float(_field(record, :maxwell_area_gate)),
        area_residual=_safe_float(_field(record, :area_residual)),
        point_requests=_safe_int(_field(_field(record, :cache_stats), :point_requests, get(cache, "point_requests", 0))),
        cache_hits=_safe_int(_field(_field(record, :cache_stats), :cache_hits, get(cache, "cache_hits", 0))),
        unique_solves=_safe_int(_field(_field(record, :cache_stats), :unique_solves, get(cache, "unique_solves", 0))),
        failed_points=failure,
        finite_and_converged=failure == 0,
    )
end

function _run_anchor(cfg, T, out_dir, telemetry)
    anchor_dir = joinpath(out_dir, "anchors", "T_$(replace(string(T), "." => "p"))")
    mkpath(anchor_dir)
    result = Models.run_production_phase_pipeline(
        :PNJL;
        T_start=T, T_end=T, dT=1.0, rho_grid=_rho_grid(), xi=cfg.xi,
        output_dir=anchor_dir,
        run_id="$(cfg.tag)_xi$(cfg.xi)_T$(T)",
        profile=:cep_hybrid_production_shadow, solver_backend=:models,
        reverse_rho=true, seed_policy=:candidates, p_num=24, t_num=8,
        thermo_quadrature_policy=:rs_reduced_adaptive,
        thermo_quadrature_rtol=1e-8, thermo_quadrature_atol=1e-10,
        thermo_quadrature_maxevals=10^7, iterations=80, compute_crossover=false,
        cep_tol=0.125, temperature_resolution_target_MeV=0.125,
        cep_max_bisect_iter=20, cep_max_refine_level_rho=0,
        # This job materializes the fixed fine pool only.  Production-parity
        # geometry is recomputed by the solver-free evaluator on selected
        # Stage-C unions, so this raw pool must not recurse to 0.0015625.
        rho_geometry_convergence=false, rho_position_tol_MeV=0.025,
        rho_density_tol=0.0025, rho_maxwell_area_tol=5e-5,
        rho_refinement_policy=:uniform_nested, rho_support_fine_step=RHO_FINE_STEP,
        rho_support_targeted_cap=12, rho_support_config=Models.RhoSupportConfig(),
        rho_hybrid_verification=Models.RhoHybridVerificationConfig(),
        work_telemetry=telemetry, memoize_uniform=true, promote_reference=false,
    )
    rows = _curve_rows(joinpath(anchor_dir, "trho_scan.csv"), cfg.xi, T, cfg.calculation_sha)
    _slice_row(result, cfg.xi, T, cfg.calculation_sha), rows
end

function _sha(path)
    bytes2hex(sha256(read(path)))
end

function _run_job(cfg)
    mkpath(cfg.output_dir)
    telemetry = Models.SolverWorkTelemetry()
    started = time_ns()
    curve_rows = NamedTuple[]
    slice_rows = NamedTuple[]
    for T in DENSITY_ANCHORS[cfg.xi]
        slice, rows = _run_anchor(cfg, T, cfg.output_dir, telemetry)
        push!(slice_rows, slice)
        append!(curve_rows, rows)
    end
    elapsed = (time_ns() - started) / 1e9
    snapshot = Models.solver_work_snapshot(telemetry)
    failed = sum(row.failed_points for row in slice_rows)
    cache = Dict(
        "point_requests" => sum(row.point_requests for row in slice_rows),
        "cache_hits" => sum(row.cache_hits for row in slice_rows),
        "unique_solves" => sum(row.unique_solves for row in slice_rows),
        "failed_points" => failed,
    )
    cache["point_request_reconciliation"] = cache["point_requests"] == cache["unique_solves"] + cache["cache_hits"]
    _write_csv(path, rows) = isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
    _write_csv(joinpath(cfg.output_dir, "fine_pool.csv"), curve_rows)
    _write_csv(joinpath(cfg.output_dir, "slice_metrics.csv"), slice_rows)
    _write_csv(joinpath(cfg.output_dir, "method_costs.csv"), [(
        xi=cfg.xi, method=METHOD, calculation_sha=cfg.calculation_sha,
        unique_solves=cache["unique_solves"], point_requests=cache["point_requests"],
        cache_hits=cache["cache_hits"], failed_points=failed,
        equilibrium_requests=snapshot.equilibrium_requests,
        fixedrho_requests=snapshot.fixedrho_requests,
        residual_calls=snapshot.nlsolve_f_calls + snapshot.postprocess_residual_calls,
        jacobian_calls=snapshot.nlsolve_g_calls, newton_iterations=snapshot.newton_iterations,
        runner_seconds=elapsed, fallback_count=snapshot.root_fallbacks,
        retry_count=snapshot.scan_retries,
    )])
    files = Dict(name => _sha(joinpath(cfg.output_dir, name)) for name in
        ("fine_pool.csv", "slice_metrics.csv", "method_costs.csv"))
    summary = Dict(
        "schema_version" => SCHEMA_VERSION, "scope" => cfg.scope, "xi" => cfg.xi,
        "method" => METHOD, "tag" => cfg.tag, "source_run_id" => cfg.source_run_id,
        "calculation_sha" => cfg.calculation_sha, "postprocess_sha" => cfg.postprocess_sha,
        "workflow_head_sha" => cfg.postprocess_sha, "rho_fine_step" => RHO_FINE_STEP,
        "rho_max" => RHO_MAX, "anchors" => DENSITY_ANCHORS[cfg.xi],
        "solver_called" => true, "finite_and_converged_final" => failed == 0,
        "cache" => cache, "telemetry" => Dict(string(field) => getproperty(snapshot, field)
            for field in propertynames(snapshot)), "runner_seconds" => elapsed,
        "files" => files,
    )
    open(joinpath(cfg.output_dir, "job_summary.json"), "w") do io
        JSON3.write(io, summary); write(io, '\n')
    end
    open(joinpath(cfg.output_dir, "manifest.json"), "w") do io
        JSON3.pretty(io, merge(summary, Dict("provenance" => Dict(
            "calculation_sha" => cfg.calculation_sha, "postprocess_sha" => cfg.postprocess_sha,
            "reference_write" => false))))
        write(io, '\n')
    end
    println(JSON3.write(summary))
end

if abspath(PROGRAM_FILE) == @__FILE__
    _run_job(_config(ARGS))
end
