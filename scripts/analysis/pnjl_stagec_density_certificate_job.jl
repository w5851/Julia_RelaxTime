#!/usr/bin/env julia

"""Run one Actions-only Stage-C density-certificate source job.

The repository checkout used for the numerical calculation is the immutable
calculation SHA.  The workflow supplies this runner from the postprocess head,
but all Models code is loaded from the calculation checkout.  This keeps the
runner/postprocess provenance separate without changing production semantics.
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
const METHODS = (:production_hybrid, :memoized_dense, :independent_oracle)
const ENDPOINT_POLICY = :three_crossing_endpoint_local_v2

@inline function _arg(args, name, default=nothing)
    index = findfirst(==(name), args)
    index === nothing && return default
    index == length(args) && throw(ArgumentError("missing value for $(name)"))
    return args[index + 1]
end

@inline function _safe_float(value)
    value === missing && return NaN
    value === nothing && return NaN
    try
        return Float64(value)
    catch
        return NaN
    end
end

@inline function _bool(value)
    value isa Bool && return value
    lowercase(strip(String(value))) in ("true", "1", "yes")
end

@inline function _field(record, name::Symbol, default=nothing)
    record === nothing && return default
    try
        return getproperty(record, name)
    catch
        return default
    end
end

function _config(args)
    xi = parse(Float64, String(_arg(args, "--xi", "0.0")))
    haskey(DENSITY_ANCHORS, xi) || throw(ArgumentError("unsupported density-feasibility xi $(xi)"))
    method = Symbol(String(_arg(args, "--method", "production_hybrid")))
    method in METHODS || throw(ArgumentError("unsupported method $(method)"))
    calculation_sha = String(_arg(args, "--calculation-sha", ""))
    occursin(r"^[0-9a-fA-F]{40}$", calculation_sha) ||
        throw(ArgumentError("calculation-sha must be a 40-character SHA"))
    output_dir = abspath(String(_arg(args, "--output-dir", joinpath(pwd(), "stagec_density_job"))))
    tag = String(_arg(args, "--tag", "cep_hybrid_stagec_density_certificate_v1"))
    postprocess_sha = String(_arg(args, "--postprocess-sha", "unknown"))
    return (; xi, method, calculation_sha, output_dir, tag, postprocess_sha)
end

function _steps(method::Symbol)
    method === :production_hybrid && return (coarse=0.05, fine=0.025, policy=:rho_support_hybrid, levels=4)
    method === :memoized_dense && return (coarse=0.0125, fine=0.00625, policy=:uniform_nested, levels=1)
    return (coarse=0.00625, fine=0.003125, policy=:uniform_nested, levels=1)
end

function _rho_grid(step::Float64)
    count = Int(round(4.0 / step))
    Float64.(collect(range(0.0, 4.0; length=count + 1)))
end

function _run_anchor(cfg, T::Float64, anchor_dir::String, steps, telemetry)
    mkpath(anchor_dir)
    Models.run_production_phase_pipeline(
        :PNJL;
        T_start=T, T_end=T, dT=1.0,
        rho_grid=_rho_grid(steps.coarse), xi=cfg.xi,
        output_dir=anchor_dir,
        run_id="$(cfg.tag)_$(cfg.method)_xi$(cfg.xi)_T$(T)",
        profile=:cep_hybrid_production_shadow,
        solver_backend=:models, reverse_rho=true, seed_policy=:candidates,
        p_num=24, t_num=8,
        thermo_quadrature_policy=:rs_reduced_adaptive,
        thermo_quadrature_rtol=1e-8, thermo_quadrature_atol=1e-10,
        thermo_quadrature_maxevals=10^7, iterations=80,
        compute_crossover=false,
        cep_tol=0.125, temperature_resolution_target_MeV=0.125,
        cep_max_bisect_iter=20,
        cep_max_refine_level_rho=steps.levels,
        rho_geometry_convergence=true,
        rho_position_tol_MeV=0.025, rho_density_tol=0.0025,
        rho_maxwell_area_tol=5e-5,
        rho_refinement_policy=steps.policy,
        rho_support_fine_step=steps.fine,
        rho_support_targeted_cap=12,
        rho_support_config=Models.RhoSupportConfig(),
        rho_hybrid_verification=Models.RhoHybridVerificationConfig(endpoint_policy=ENDPOINT_POLICY),
        work_telemetry=telemetry,
        memoize_uniform=steps.policy === :uniform_nested,
        promote_reference=false,
    )
end

function _anchor_record(diag, T::Float64)
    records = get(diag, "sweep_records", NamedTuple[])
    for record in records
        candidate = _safe_float(_field(record, :T_MeV, NaN))
        isfinite(candidate) && isapprox(candidate, T; atol=1e-8, rtol=0.0) && return record
    end
    nothing
end

function _curve_rows(path::String, cfg, steps)
    isfile(path) || return NamedTuple[]
    rows = NamedTuple[]
    seen = Set{Tuple{Float64, Float64, Float64}}()
    for row in CSV.File(path)
        T = _safe_float(_field(row, :T_MeV, NaN))
        rho = _safe_float(_field(row, :rho, NaN))
        xi = _safe_float(_field(row, :xi, NaN))
        isfinite(T) && isfinite(rho) && isfinite(xi) || continue
        key = (T, xi, rho)
        key in seen && continue
        push!(seen, key)
        on_coarse = isapprox(rho / steps.coarse, round(rho / steps.coarse); atol=1e-8, rtol=0.0)
        on_fine = isapprox(rho / steps.fine, round(rho / steps.fine); atol=1e-8, rtol=0.0)
        role, level = if cfg.method === :production_hybrid
            on_coarse ? ("coarse_grid", 0) : on_fine ? ("fine_grid", 1) :
            ("stage_c_guard", 4)
        elseif on_coarse
            ("coarse_grid", 0)
        elseif on_fine
            ("fine_grid", 1)
        else
            ("targeted", 2)
        end
        converged = _bool(_field(row, :converged, false))
        residual = _safe_float(_field(row, :residual_norm, NaN))
        push!(rows, (
            xi=xi, method=String(cfg.method), calculation_sha=cfg.calculation_sha,
            T_MeV=T, rho_level=level, rho=rho,
            muq_MeV=_safe_float(_field(row, :mu_avg_MeV, NaN)),
            pressure_fm4=_safe_float(_field(row, :pressure_fm4, NaN)),
            residual_norm=residual, iterations=Int(round(_safe_float(_field(row, :iterations, 0)))),
            converged=converged, finite=converged && isfinite(residual), sampling_role=role,
            cache_reused=false, solver_message=String(_field(row, :message, "")),
        ))
    end
    sort!(rows; by=row -> (row.T_MeV, row.rho_level, row.rho))
    rows
end

function _slice_row(cfg, T, result)
    diag = result.diagnostics
    record = _anchor_record(diag, T)
    cache = get(diag, "rho_support_cache", Dict{String, Any}())
    cache_dict = cache isa AbstractDict ? cache : Dict{String, Any}()
    failure = Int(_field(record, :stats_failure, get(diag, "scan_failure", 0)))
    endpoint_left = _field(record, :hybrid_endpoint_left_bracket, nothing)
    endpoint_right = _field(record, :hybrid_endpoint_right_bracket, nothing)
    endpoint_low = _safe_float(_field(endpoint_left, :low, NaN))
    endpoint_high = _safe_float(_field(endpoint_left, :high, NaN))
    endpoint_right_low = _safe_float(_field(endpoint_right, :low, NaN))
    endpoint_right_high = _safe_float(_field(endpoint_right, :high, NaN))
    status = String(_field(record, :status, result.cep.result_status))
    return (
        xi=cfg.xi, method=String(cfg.method), calculation_sha=cfg.calculation_sha, T_MeV=T,
        stage_a_status=String(_field(record, :hybrid_stage_a_status, _field(record, :cascade_status, :not_run))),
        stage_b_status=String(_field(record, :hybrid_stage_b_status, :not_run)),
        stage_c_status=String(_field(record, :hybrid_stage_c_status, :not_run)),
        stage_used=String(_field(record, :stage_used, :uniform_nested)),
        result_status=status, raw_status=String(_field(record, :raw_status, status)),
        geometry_converged=Bool(_field(record, :geometry_converged, false)),
        position_error_MeV=_safe_float(_field(record, :position_error_MeV, NaN)),
        density_error=_safe_float(_field(record, :density_error, NaN)),
        maxwell_area_gate=_safe_float(_field(record, :maxwell_area_gate, NaN)),
        area_residual=_safe_float(_field(record, :area_residual, NaN)),
        rho_hadron=_safe_float(_field(record, :rho_hadron, NaN)),
        rho_quark=_safe_float(_field(record, :rho_quark, NaN)),
        mu_transition_MeV=_safe_float(_field(record, :mu_transition_MeV, _field(record, :mu_transition, NaN))),
        rho_spinodal_hadron=_safe_float(_field(record, :rho_spinodal_hadron, NaN)),
        rho_spinodal_quark=_safe_float(_field(record, :rho_spinodal_quark, NaN)),
        support_low=_safe_float(_field(record, :hybrid_support_low, NaN)),
        support_high=_safe_float(_field(record, :hybrid_support_high, NaN)),
        endpoint_left_bracket_low=endpoint_low, endpoint_left_bracket_high=endpoint_high,
        endpoint_right_bracket_low=endpoint_right_low, endpoint_right_bracket_high=endpoint_right_high,
        targeted_additions=Int(_field(record, :cascade_targeted_count, get(cache_dict, "targeted_additions", 0))),
        stage_c_targeted_additions=Int(_field(record, :hybrid_targeted_point_count, 0)),
        unique_solves=Int(_field(_field(record, :cache_stats, nothing), :unique_solves, get(cache_dict, "unique_solves", 0))),
        point_requests=Int(_field(_field(record, :cache_stats, nothing), :point_requests, get(cache_dict, "point_requests", 0))),
        cache_hits=Int(_field(_field(record, :cache_stats, nothing), :cache_hits, get(cache_dict, "cache_hits", 0))),
        finite_and_converged=failure == 0, solver_failure_count=failure,
    )
end

function _write_csv(path, rows)
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _run_job(cfg)
    steps = _steps(cfg.method)
    telemetry = Models.SolverWorkTelemetry()
    started = time_ns()
    curve_rows = NamedTuple[]
    slice_rows = NamedTuple[]
    accuracy_rows = NamedTuple[]
    for T in DENSITY_ANCHORS[cfg.xi]
        anchor_dir = joinpath(cfg.output_dir, "anchors", "T_$(replace(string(T), "." => "p"))")
        result = _run_anchor(cfg, T, anchor_dir, steps, telemetry)
        append!(curve_rows, _curve_rows(joinpath(anchor_dir, "trho_scan.csv"), cfg, steps))
        push!(slice_rows, _slice_row(cfg, T, result))
        cep = result.cep
        push!(accuracy_rows, (
            xi=cfg.xi, method=String(cfg.method), calculation_sha=cfg.calculation_sha,
            anchor_T_MeV=T, result_status=String(cep.result_status), found=cep.found,
            T_last_first_order_MeV=cep.T_last_first_order_MeV,
            T_first_monotone_MeV=cep.T_first_monotone_MeV,
            ambiguity_width_T_MeV=cep.ambiguity_width_T_MeV,
            bracket_width_T_MeV=cep.ambiguity_width_T_MeV,
            temperature_resolution_target_MeV=cep.temperature_resolution_target_MeV,
        ))
    end
    elapsed = (time_ns() - started) / 1e9
    snapshot = Models.solver_work_snapshot(telemetry)
    cache = Dict(
        "point_requests" => sum(Int(row.point_requests) for row in slice_rows),
        "cache_hits" => sum(Int(row.cache_hits) for row in slice_rows),
        "unique_solves" => sum(Int(row.unique_solves) for row in slice_rows),
        "targeted_additions" => sum(Int(row.targeted_additions) for row in slice_rows),
        "failed_points" => sum(Int(row.solver_failure_count) for row in slice_rows),
    )
    cache["point_request_reconciliation"] = cache["point_requests"] == cache["unique_solves"] + cache["cache_hits"]
    mkpath(cfg.output_dir)
    _write_csv(joinpath(cfg.output_dir, "curve_points.csv"), curve_rows)
    _write_csv(joinpath(cfg.output_dir, "slice_metrics.csv"), slice_rows)
    _write_csv(joinpath(cfg.output_dir, "cep_accuracy.csv"), accuracy_rows)
    _write_csv(joinpath(cfg.output_dir, "method_costs.csv"), [(
        xi=cfg.xi, method=String(cfg.method), calculation_sha=cfg.calculation_sha,
        unique_solves=cache["unique_solves"], equilibrium_requests=snapshot.equilibrium_requests,
        fixedrho_requests=snapshot.fixedrho_requests, requested_point_calls=cache["point_requests"],
        cache_hits=cache["cache_hits"], targeted_additions=cache["targeted_additions"],
        residual_calls=snapshot.nlsolve_f_calls + snapshot.postprocess_residual_calls,
        jacobian_calls=snapshot.nlsolve_g_calls, newton_iterations=snapshot.newton_iterations,
        runner_seconds=elapsed, fallback_count=snapshot.root_fallbacks, retry_count=snapshot.scan_retries,
    )])
    hashes = Dict(name => bytes2hex(sha256(read(joinpath(cfg.output_dir, name)))) for name in ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv"))
    summary = Dict(
        "schema_version" => "cep_hybrid_stagec_density_certificate_job_v1",
        "xi" => cfg.xi, "method" => String(cfg.method), "tag" => cfg.tag,
        "calculation_sha" => cfg.calculation_sha, "postprocess_sha" => cfg.postprocess_sha,
        "workflow_head_sha" => cfg.postprocess_sha,
        "scope" => "density_feasibility", "anchors" => DENSITY_ANCHORS[cfg.xi],
        "parameters" => Dict("stage_c_local_step" => STAGE_C_STEP, "rho_support_targeted_cap" => 12,
            "rho_refinement_policy" => String(steps.policy),
            "rho_point_ranking_version" => steps.policy === :rho_support_hybrid ? "stage_b_features_v1" : "uniform_nested_reference",
            "rho_hybrid_endpoint_policy" => String(ENDPOINT_POLICY)),
        "telemetry" => Dict(string(field) => getproperty(snapshot, field) for field in propertynames(snapshot)),
        "cache" => cache, "runner_seconds" => elapsed, "finite_and_converged_final" => all(row -> row.finite_and_converged, slice_rows),
        "curve_file_sha256" => hashes, "provenance" => Dict("calculation_sha" => cfg.calculation_sha,
            "postprocess_sha" => cfg.postprocess_sha, "reference_write" => false),
    )
    open(joinpath(cfg.output_dir, "job_summary.json"), "w") do io
        JSON3.write(io, summary); write(io, '\n')
    end
    open(joinpath(cfg.output_dir, "manifest.json"), "w") do io
        JSON3.pretty(io, Dict("schema_version" => summary["schema_version"], "calculation_sha" => cfg.calculation_sha,
            "postprocess_sha" => cfg.postprocess_sha, "workflow_head_sha" => cfg.postprocess_sha,
            "method" => String(cfg.method), "xi" => cfg.xi,
            "scope" => "density_feasibility", "files" => hashes)); write(io, '\n')
    end
    println(JSON3.write(summary))
end

if abspath(PROGRAM_FILE) == @__FILE__
    _run_job(_config(ARGS))
end
