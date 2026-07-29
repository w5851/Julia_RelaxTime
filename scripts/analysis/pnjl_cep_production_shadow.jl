#!/usr/bin/env julia

"""Run one fixed-anchor PNJL CEP production-shadow job.

This runner is intentionally opt-in and diagnostic-only.  It exercises the
production phase pipeline for the cascade method and request-scoped memoized
uniform sessions for the two comparison methods.  It never writes the
reference tree or launches a full phase-reference/transport workflow.
"""

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using CSV
using Dates
using JSON3
using SHA

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
end
using Main.Models

const SHADOW_XI = (-0.5, 0.0, 0.5)
const SHADOW_ANCHORS = Dict(
    -0.5 => [5.0, 20.0, 60.0, 100.0, 130.0, 147.0947265625, 147.2197265625, 160.0],
    0.0 => [5.0, 20.0, 60.0, 100.0, 120.0, 130.9619140625, 131.0869140625, 145.0],
    0.5 => [5.0, 20.0, 60.0, 90.0, 100.0, 106.9599609375, 107.0849609375, 120.0],
)
const SHADOW_METHODS = (:production_cascade, :memoized_dense, :independent_oracle)
const HBARC_SHADOW = 197.3269804

@inline function _arg(values, name, default=nothing)
    index = findfirst(==(name), values)
    index === nothing && return default
    index == length(values) && throw(ArgumentError("missing value for $(name)"))
    return values[index + 1]
end

function _config(args)
    method = Symbol(_arg(args, "--method", "production_cascade"))
    method in SHADOW_METHODS || throw(ArgumentError("unsupported shadow method $(method)"))
    xi = parse(Float64, _arg(args, "--xi", "0.0"))
    any(isapprox(xi, value; atol=1e-9, rtol=0.0) for value in SHADOW_XI) ||
        throw(ArgumentError("shadow xi must be one of -0.5, 0.0, 0.5"))
    calculation_sha = _arg(args, "--calculation-sha", "")
    occursin(r"^[0-9a-fA-F]{40}$", calculation_sha) ||
        throw(ArgumentError("calculation-sha must be an immutable 40-character SHA"))
    output_dir = abspath(_arg(args, "--output-dir", joinpath(pwd(), "shadow_artifact")))
    tag = String(_arg(args, "--tag", "cep_cascade_production_shadow_v1"))
    return (; method, xi, calculation_sha, output_dir, tag)
end

@inline function _method_steps(method::Symbol)
    method === :production_cascade && return (coarse=0.05, fine=0.025, policy=:rho_support_cascade)
    method === :memoized_dense && return (coarse=0.0125, fine=0.00625, policy=:uniform_nested)
    return (coarse=0.00625, fine=0.003125, policy=:uniform_nested)
end

function _rho_grid(step::Float64)
    count = Int(round(4.0 / step))
    return Float64.(collect(range(0.0, 4.0; length=count + 1)))
end

@inline function _bool(value)
    value isa Bool && return value
    lowercase(strip(String(value))) in ("true", "1", "yes")
end

function _safe_float(value)
    value === missing && return NaN
    value === nothing && return NaN
    try
        return Float64(value)
    catch
        return NaN
    end
end

function _curve_rows(path::String, cfg, steps)
    isfile(path) || return NamedTuple[]
    rows = NamedTuple[]
    seen = Set{Tuple{Float64, Float64, Float64} }
    for row in CSV.File(path)
        T = _safe_float(getproperty(row, :T_MeV))
        rho = _safe_float(getproperty(row, :rho))
        xi = _safe_float(getproperty(row, :xi))
        isfinite(T) && isfinite(rho) && isfinite(xi) || continue
        key = (T, xi, rho)
        key in seen && continue
        push!(seen, key)
        on_fine = isapprox(rho / steps.fine, round(rho / steps.fine); atol=1e-8, rtol=0.0)
        on_coarse = isapprox(rho / steps.coarse, round(rho / steps.coarse); atol=1e-8, rtol=0.0)
        role = on_coarse || on_fine ? "grid" : "targeted"
        push!(rows, (
            xi=xi,
            method=String(cfg.method),
            calculation_sha=cfg.calculation_sha,
            T_MeV=T,
            rho_level=on_coarse ? 0 : 1,
            rho=rho,
            muq_MeV=_safe_float(getproperty(row, :mu_avg_MeV)),
            pressure_fm4=_safe_float(getproperty(row, :pressure_fm4)),
            residual_norm=_safe_float(getproperty(row, :residual_norm)),
            iterations=Int(round(_safe_float(getproperty(row, :iterations)))),
            converged=_bool(getproperty(row, :converged)),
            finite=_bool(getproperty(row, :converged)) && isfinite(_safe_float(getproperty(row, :residual_norm))),
            sampling_role=role,
            cache_reused=false,
            solver_message=String(getproperty(row, :message)),
        ))
    end
    sort!(rows; by=row -> (row.T_MeV, row.rho_level, row.rho))
    return rows
end

function _write_csv(path, rows)
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _run_anchor(cfg, T::Float64, anchor_dir::String, steps, telemetry)
    mkpath(anchor_dir)
    memoize_uniform = cfg.method !== :production_cascade
    result = Models.run_production_phase_pipeline(
        :PNJL;
        T_start=T,
        T_end=T,
        dT=1.0,
        rho_grid=_rho_grid(steps.coarse),
        xi=cfg.xi,
        output_dir=anchor_dir,
        run_id="$(cfg.tag)_$(cfg.method)_xi$(cfg.xi)_T$(T)",
        profile=:cep_production_shadow,
        solver_backend=:models,
        reverse_rho=true,
        seed_policy=:candidates,
        p_num=24,
        t_num=8,
        thermo_quadrature_policy=:rs_reduced_adaptive,
        thermo_quadrature_rtol=1e-8,
        thermo_quadrature_atol=1e-10,
        thermo_quadrature_maxevals=10^7,
        iterations=80,
        compute_crossover=false,
        cep_tol=0.125,
        temperature_resolution_target_MeV=0.125,
        cep_max_bisect_iter=20,
        cep_max_refine_level_rho=1,
        rho_geometry_convergence=true,
        rho_position_tol_MeV=0.025,
        rho_density_tol=0.0025,
        rho_maxwell_area_tol=5e-5,
        rho_refinement_policy=steps.policy,
        rho_support_fine_step=steps.fine,
        rho_support_targeted_cap=12,
        rho_support_config=Models.RhoSupportConfig(),
        work_telemetry=telemetry,
        memoize_uniform=memoize_uniform,
        promote_reference=false,
    )
    return result
end

function _anchor_record(diag, T::Float64)
    records = get(diag, "sweep_records", NamedTuple[])
    for record in records
        candidate = try
            Float64(getproperty(record, :T_MeV))
        catch
            NaN
        end
        isfinite(candidate) && isapprox(candidate, T; atol=1e-8, rtol=0.0) && return record
    end
    return nothing
end

function _write_job(cfg)
    steps = _method_steps(cfg.method)
    started = time_ns()
    telemetry = Models.SolverWorkTelemetry()
    curve_rows = NamedTuple[]
    slice_rows = NamedTuple[]
    accuracy_rows = NamedTuple[]
    anchor_dirs = String[]
    results = Models.PhasePipelineResult[]

    for T in SHADOW_ANCHORS[cfg.xi]
        anchor_dir = joinpath(cfg.output_dir, "anchors", "T_$(replace(string(T), "." => "p"))")
        push!(anchor_dirs, anchor_dir)
        result = _run_anchor(cfg, T, anchor_dir, steps, telemetry)
        push!(results, result)
        append!(curve_rows, _curve_rows(joinpath(anchor_dir, "trho_scan.csv"), cfg, steps))
        diag = result.diagnostics
        statuses = get(diag, "sweep_statuses", String[])
        reasons = get(diag, "sweep_reasons", String[])
        record = _anchor_record(diag, T)
        status = record === nothing ? (isempty(statuses) ? String(result.cep.result_status) : String(first(statuses))) :
            String(get(record, :status, result.cep.result_status))
        reason = record === nothing ? (isempty(reasons) ? "" : String(first(reasons))) :
            String(get(record, :reason, ""))
        cache = get(diag, "rho_support_cache", nothing)
        cache_dict = cache === nothing ? Dict{String, Any}() : cache
        geometry_converged = record === nothing ? get(diag, "rho_unconverged_count", 0) == 0 :
            Bool(get(record, :geometry_converged, false))
        cache_stats = record === nothing ? nothing : get(record, :cache_stats, nothing)
        targeted_fallback = cache_stats === nothing ? 0 : Int(get(cache_stats, :targeted_additions, 0))
        targeted_additions = record === nothing ? 0 : Int(get(record, :cascade_targeted_count, targeted_fallback))
        unique_solves = cache_stats === nothing ? Int(get(cache_dict, "unique_solves", get(diag, "scan_total", 0))) :
            Int(get(cache_stats, :unique_solves, 0))
        point_requests = cache_stats === nothing ? Int(get(cache_dict, "point_requests", get(diag, "scan_total", 0))) :
            Int(get(cache_stats, :point_requests, 0))
        cache_hits_for_slice = cache_stats === nothing ? Int(get(cache_dict, "cache_hits", 0)) :
            Int(get(cache_stats, :cache_hits, 0))
        solver_failure_count = record === nothing ? Int(get(diag, "scan_failure", 0)) :
            Int(get(record, :stats_failure, 0))
        push!(slice_rows, (
            xi=cfg.xi,
            method=String(cfg.method),
            calculation_sha=cfg.calculation_sha,
            T_MeV=T,
            coarse_status=record === nothing ? status : String(get(record, :coarse_status, status)),
            fine_status=record === nothing ? status : String(get(record, :fine_status, status)),
            result_status=status,
            coarse_reason=record === nothing ? reason : String(get(record, :coarse_reason, reason)),
            fine_reason=record === nothing ? reason : String(get(record, :fine_reason, reason)),
            cascade_status=record === nothing ? "not_run" : String(get(record, :cascade_status, "not_run")),
            geometry_converged=geometry_converged,
            position_error_MeV=record === nothing ? NaN : Float64(get(record, :position_error_MeV, NaN)),
            density_error=record === nothing ? NaN : Float64(get(record, :density_error, NaN)),
            maxwell_area_gate=record === nothing ? NaN : Float64(get(record, :maxwell_area_gate, NaN)),
            area_residual=record === nothing ? NaN : Float64(get(record, :area_residual, NaN)),
            rho_hadron=record === nothing ? NaN : Float64(get(record, :rho_hadron, NaN)),
            rho_quark=record === nothing ? NaN : Float64(get(record, :rho_quark, NaN)),
            mu_spinodal_hadron_MeV=record === nothing ? NaN : Float64(get(record, :mu_spinodal_hadron_MeV, NaN)),
            mu_spinodal_quark_MeV=record === nothing ? NaN : Float64(get(record, :mu_spinodal_quark_MeV, NaN)),
            rho_spinodal_hadron=record === nothing ? NaN : Float64(get(record, :rho_spinodal_hadron, NaN)),
            rho_spinodal_quark=record === nothing ? NaN : Float64(get(record, :rho_spinodal_quark, NaN)),
            targeted_additions=targeted_additions,
            unique_solves=unique_solves,
            point_requests=point_requests,
            cache_hits=cache_hits_for_slice,
            finite_and_converged=Bool(solver_failure_count == 0),
            solver_failure_count=solver_failure_count,
        ))
        cep = result.cep
        push!(accuracy_rows, (
            xi=cfg.xi,
            method=String(cfg.method),
            calculation_sha=cfg.calculation_sha,
            anchor_T_MeV=T,
            result_status=String(cep.result_status),
            found=cep.found,
            T_last_first_order_MeV=cep.T_last_first_order_MeV,
            muq_last_first_order_MeV=cep.mu_last_first_order_MeV,
            T_first_monotone_MeV=cep.T_first_monotone_MeV,
            ambiguity_width_T_MeV=cep.ambiguity_width_T_MeV,
            temperature_resolution_target_MeV=cep.temperature_resolution_target_MeV,
        ))
    end

    snapshot = Models.solver_work_snapshot(telemetry)
    elapsed = (time_ns() - started) / 1e9
    cache = get(first(results).diagnostics, "rho_support_cache", nothing)
    cache_dict = cache === nothing ? Dict{String, Any}() : cache
    costs = [(
        xi=cfg.xi,
        method=String(cfg.method),
        calculation_sha=cfg.calculation_sha,
        equilibrium_requests=snapshot.equilibrium_requests,
        fixedrho_requests=snapshot.fixedrho_requests,
        unique_solves=Int(get(cache_dict, "unique_solves", sum(row.unique_solves for row in slice_rows))),
        requested_point_calls=Int(get(cache_dict, "point_requests", sum(row.point_requests for row in slice_rows))),
        uncached_equivalent_requests=Int(get(cache_dict, "point_requests", sum(row.point_requests for row in slice_rows))),
        cache_hits=Int(get(cache_dict, "cache_hits", sum(row.cache_hits for row in slice_rows))),
        targeted_additions=Int(get(cache_dict, "targeted_additions", sum(row.targeted_additions for row in slice_rows))),
        residual_calls=snapshot.nlsolve_f_calls + snapshot.postprocess_residual_calls,
        jacobian_calls=snapshot.nlsolve_g_calls,
        newton_iterations=snapshot.newton_iterations,
        trust_region_iterations=snapshot.trust_region_iterations,
        nonconverged_attempts=snapshot.nonconverged_attempts,
        fallback_count=snapshot.root_fallbacks,
        governed_rescue_count=snapshot.governed_rescues,
        retry_count=snapshot.scan_retries,
        exception_count=snapshot.exceptions,
        runner_seconds=elapsed,
    )]
    _write_csv(joinpath(cfg.output_dir, "curve_points.csv"), curve_rows)
    _write_csv(joinpath(cfg.output_dir, "slice_metrics.csv"), slice_rows)
    _write_csv(joinpath(cfg.output_dir, "method_costs.csv"), costs)
    _write_csv(joinpath(cfg.output_dir, "cep_accuracy.csv"), accuracy_rows)

    hashes = Dict{String, String}()
    for path in filter(isfile, [joinpath(cfg.output_dir, name) for name in ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv")])
        hashes[basename(path)] = bytes2hex(sha256(read(path)))
    end
    summary = Dict(
        "schema_version" => "cep_cascade_production_shadow_v1",
        "xi" => cfg.xi,
        "method" => String(cfg.method),
        "tag" => cfg.tag,
        "calculation_sha" => cfg.calculation_sha,
        "workflow_head_sha" => cfg.calculation_sha,
        "anchors" => SHADOW_ANCHORS[cfg.xi],
        "parameters" => Dict("p_num" => 24, "t_num" => 8,
            "thermo_quadrature_policy" => "rs_reduced_adaptive",
            "thermo_quadrature_rtol" => 1e-8, "thermo_quadrature_atol" => 1e-10,
            "thermo_quadrature_maxevals" => 10^7, "iterations" => 80,
            "rho_coarse_step" => steps.coarse, "rho_fine_step" => steps.fine,
            "rho_support_targeted_cap" => cfg.method === :production_cascade ? 12 : 0),
        "telemetry" => Dict(string(field) => getproperty(snapshot, field) for field in propertynames(snapshot)),
        "cache" => cache_dict,
        "runner_seconds" => elapsed,
        "finite_and_converged_final" => all(row -> row.finite_and_converged, slice_rows),
        "curve_file_sha256" => hashes,
        "provenance" => Dict("calculation_sha" => cfg.calculation_sha,
            "reference_write" => false, "anchor_run_count" => length(anchor_dirs)),
    )
    open(joinpath(cfg.output_dir, "job_summary.json"), "w") do io
        JSON3.write(io, summary)
        write(io, '\n')
    end
    open(joinpath(cfg.output_dir, "manifest.json"), "w") do io
        JSON3.pretty(io, Dict("schema_version" => "cep_cascade_production_shadow_v1",
            "calculation_sha" => cfg.calculation_sha, "files" => hashes,
            "anchors" => SHADOW_ANCHORS[cfg.xi], "method" => String(cfg.method)))
    end
    println(JSON3.write(summary))
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    _write_job(_config(ARGS))
end
