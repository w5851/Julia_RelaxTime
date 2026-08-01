#!/usr/bin/env julia

"""Run one independent high-resolution PNJL deep-oracle anchor.

This is a diagnostic-only runner.  It keeps the production calculation SHA
separate from the workflow/post-processing SHA and never writes reference
artifacts.  The runner reuses the existing production-shadow anchor helpers,
but fixes the independent rho grid to 0.003125 -> 0.0015625.
"""

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using CSV
using JSON3
using SHA

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
end
using Main.Models

# Reuse the production-shadow implementation of the PhaseCore/Maxwell path.
# Its main entrypoint is guarded by PROGRAM_FILE, so including it does not
# launch the normal nine-anchor job.
include(joinpath(@__DIR__, "pnjl_cep_production_shadow.jl"))

const DEEP_COARSE_STEP = 0.003125
const DEEP_FINE_STEP = 0.0015625
const DEEP_METHOD = :independent_oracle

@inline function _deep_arg(values, name, default=nothing)
    index = findfirst(==(name), values)
    index === nothing && return default
    index == length(values) && throw(ArgumentError("missing value for $(name)"))
    return values[index + 1]
end

function _deep_config(args)
    xi = parse(Float64, _deep_arg(args, "--xi", "0.0"))
    T = parse(Float64, _deep_arg(args, "--temperature", "20.0"))
    calculation_sha = String(_deep_arg(args, "--calculation-sha", ""))
    workflow_head_sha = String(_deep_arg(args, "--workflow-head-sha", ""))
    occursin(r"^[0-9a-fA-F]{40}$", calculation_sha) ||
        throw(ArgumentError("calculation-sha must be an immutable 40-character SHA"))
    occursin(r"^[0-9a-fA-F]{40}$", workflow_head_sha) ||
        throw(ArgumentError("workflow-head-sha must be an immutable 40-character SHA"))
    isfinite(xi) && -0.5 <= xi <= 0.5 || throw(ArgumentError("xi must be finite and in [-0.5, 0.5]"))
    isfinite(T) && T > 0 || throw(ArgumentError("temperature must be finite and positive"))
    output_dir = abspath(String(_deep_arg(args, "--output-dir", joinpath(pwd(), "deep_oracle_artifact"))))
    tag = String(_deep_arg(args, "--tag", "cep_deep_oracle_v1"))
    return (; xi, T, calculation_sha, workflow_head_sha, output_dir, tag)
end

function _deep_write_csv(path, rows)
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _deep_record_value(record, key::Symbol, default)
    record === nothing && return default
    try
        return get(record, key, default)
    catch
        try
            return getproperty(record, key)
        catch
            return default
        end
    end
end

function _deep_finite(value)
    try
        return isfinite(Float64(value)) ? Float64(value) : nothing
    catch
        return nothing
    end
end

function _deep_write_job(cfg)
    mkpath(cfg.output_dir)
    started = time_ns()
    telemetry = Models.SolverWorkTelemetry()
    steps = (coarse=DEEP_COARSE_STEP, fine=DEEP_FINE_STEP, policy=:uniform_nested)
    anchor_dir = joinpath(cfg.output_dir, "anchor", "T_$(replace(string(cfg.T), "." => "p"))")
    result = _run_anchor((method=DEEP_METHOD, xi=cfg.xi,
            calculation_sha=cfg.calculation_sha, output_dir=cfg.output_dir, tag=cfg.tag),
        cfg.T, anchor_dir, steps, telemetry)

    curve_rows = _curve_rows(joinpath(anchor_dir, "trho_scan.csv"),
        (method=DEEP_METHOD, xi=cfg.xi, calculation_sha=cfg.calculation_sha), steps)
    diag = result.diagnostics
    record = _anchor_record(diag, cfg.T)
    cache = get(diag, "rho_support_cache", Dict{String, Any}())
    cache isa AbstractDict || (cache = Dict{String, Any}())
    status = String(_deep_record_value(record, :status, result.cep.result_status))
    geometry_converged = Bool(_deep_record_value(record, :geometry_converged,
        get(diag, "rho_unconverged_count", 0) == 0))
    cache_stats = _deep_record_value(record, :cache_stats, nothing)
    unique_solves = Int(_deep_record_value(cache_stats, :unique_solves,
        get(cache, "unique_solves", get(diag, "scan_total", 0))))
    point_requests = Int(_deep_record_value(cache_stats, :point_requests,
        get(cache, "point_requests", get(diag, "scan_total", 0))))
    cache_hits = Int(_deep_record_value(cache_stats, :cache_hits,
        get(cache, "cache_hits", 0)))
    solver_failures = Int(_deep_record_value(record, :stats_failure,
        get(diag, "scan_failure", 0)))
    targeted = Int(_deep_record_value(record, :cascade_targeted_count, 0))
    slice_rows = [(
        xi=cfg.xi,
        method=String(DEEP_METHOD),
        calculation_sha=cfg.calculation_sha,
        workflow_head_sha=cfg.workflow_head_sha,
        T_MeV=cfg.T,
        coarse_step=DEEP_COARSE_STEP,
        fine_step=DEEP_FINE_STEP,
        coarse_status=String(_deep_record_value(record, :coarse_status, status)),
        fine_status=String(_deep_record_value(record, :fine_status, status)),
        result_status=status,
        coarse_reason=String(_deep_record_value(record, :coarse_reason, "")),
        fine_reason=String(_deep_record_value(record, :fine_reason, "")),
        geometry_converged=geometry_converged,
        position_error_MeV=Float64(_deep_record_value(record, :position_error_MeV, NaN)),
        density_error=Float64(_deep_record_value(record, :density_error, NaN)),
        maxwell_area_gate=Float64(_deep_record_value(record, :maxwell_area_gate, NaN)),
        area_residual=Float64(_deep_record_value(record, :area_residual, NaN)),
        rho_hadron=Float64(_deep_record_value(record, :rho_hadron, NaN)),
        rho_quark=Float64(_deep_record_value(record, :rho_quark, NaN)),
        mu_spinodal_hadron_MeV=Float64(_deep_record_value(record, :mu_spinodal_hadron_MeV, NaN)),
        mu_spinodal_quark_MeV=Float64(_deep_record_value(record, :mu_spinodal_quark_MeV, NaN)),
        rho_spinodal_hadron=Float64(_deep_record_value(record, :rho_spinodal_hadron, NaN)),
        rho_spinodal_quark=Float64(_deep_record_value(record, :rho_spinodal_quark, NaN)),
        targeted_additions=targeted,
        unique_solves=unique_solves,
        point_requests=point_requests,
        cache_hits=cache_hits,
        finite_and_converged=solver_failures == 0,
        solver_failure_count=solver_failures,
    )]

    snapshot = Models.solver_work_snapshot(telemetry)
    elapsed = (time_ns() - started) / 1e9
    costs = [(
        xi=cfg.xi,
        method=String(DEEP_METHOD),
        calculation_sha=cfg.calculation_sha,
        workflow_head_sha=cfg.workflow_head_sha,
        T_MeV=cfg.T,
        equilibrium_requests=snapshot.equilibrium_requests,
        fixedrho_requests=snapshot.fixedrho_requests,
        unique_solves=unique_solves,
        requested_point_calls=point_requests,
        cache_hits=cache_hits,
        targeted_additions=targeted,
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

    cep = result.cep
    accuracy = [(
        xi=cfg.xi,
        method=String(DEEP_METHOD),
        calculation_sha=cfg.calculation_sha,
        workflow_head_sha=cfg.workflow_head_sha,
        T_MeV=cfg.T,
        result_status=String(cep.result_status),
        found=cep.found,
        T_last_first_order_MeV=_deep_finite(cep.T_last_first_order_MeV),
        muq_last_first_order_MeV=_deep_finite(cep.muq_last_first_order_MeV),
        T_first_monotone_MeV=_deep_finite(cep.T_first_monotone_MeV),
        ambiguity_width_T_MeV=_deep_finite(cep.ambiguity_width_T_MeV),
        temperature_resolution_target_MeV=cep.temperature_resolution_target_MeV,
    )]

    files = Dict{String, String}()
    _deep_write_csv(joinpath(cfg.output_dir, "curve_points.csv"), curve_rows)
    _deep_write_csv(joinpath(cfg.output_dir, "slice_metrics.csv"), slice_rows)
    _deep_write_csv(joinpath(cfg.output_dir, "method_costs.csv"), costs)
    _deep_write_csv(joinpath(cfg.output_dir, "cep_accuracy.csv"), accuracy)
    for name in ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv")
        path = joinpath(cfg.output_dir, name)
        files[name] = bytes2hex(sha256(read(path)))
    end

    summary = Dict{String, Any}(
        "schema_version" => "cep_deep_oracle_v1",
        "xi" => cfg.xi,
        "temperature_MeV" => cfg.T,
        "method" => String(DEEP_METHOD),
        "tag" => cfg.tag,
        "calculation_sha" => cfg.calculation_sha,
        "workflow_head_sha" => cfg.workflow_head_sha,
        "rho_coarse_step" => DEEP_COARSE_STEP,
        "rho_fine_step" => DEEP_FINE_STEP,
        "result_status" => String(cep.result_status),
        "geometry_converged" => geometry_converged,
        "finite_and_converged_final" => solver_failures == 0,
        "unique_solves" => unique_solves,
        "point_requests" => point_requests,
        "cache_hits" => cache_hits,
        "runner_seconds" => elapsed,
        "telemetry" => Dict(string(field) => getproperty(snapshot, field) for field in propertynames(snapshot)),
        "curve_file_sha256" => files,
        "provenance" => Dict(
            "calculation_sha" => cfg.calculation_sha,
            "workflow_head_sha" => cfg.workflow_head_sha,
            "reference_write" => false,
            "solver_session" => "independent_single_anchor",
        ),
    )
    open(joinpath(cfg.output_dir, "job_summary.json"), "w") do io
        JSON3.write(io, summary)
        write(io, '\n')
    end
    manifest = Dict(
        "schema_version" => "cep_deep_oracle_v1",
        "calculation_sha" => cfg.calculation_sha,
        "workflow_head_sha" => cfg.workflow_head_sha,
        "xi" => cfg.xi,
        "temperature_MeV" => cfg.T,
        "rho_coarse_step" => DEEP_COARSE_STEP,
        "rho_fine_step" => DEEP_FINE_STEP,
        "files" => files,
    )
    open(joinpath(cfg.output_dir, "manifest.json"), "w") do io
        JSON3.pretty(io, manifest)
        write(io, '\n')
    end
    println(JSON3.write(summary))
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    _deep_write_job(_deep_config(ARGS))
end
