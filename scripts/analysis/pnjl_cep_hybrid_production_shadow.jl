#!/usr/bin/env julia

"""Run one fixed-anchor PNJL CEP hybrid-production shadow job.

The job is deliberately diagnostic-only.  It exercises the opt-in hybrid
rho policy (or one of its independent comparison methods), writes immutable
job-local evidence, and never touches the reference tree.
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

const HYBRID_XI = (-0.5, 0.0, 0.5)
const HYBRID_ANCHORS = Dict(
    -0.5 => [5.0, 20.0, 60.0, 100.0, 130.0, 147.0947265625, 147.2197265625, 160.0],
    0.0 => [5.0, 20.0, 60.0, 100.0, 120.0, 130.9619140625, 131.0869140625, 145.0],
    0.5 => [5.0, 20.0, 60.0, 90.0, 100.0, 106.9599609375, 107.0849609375, 120.0],
)
const HYBRID_METHODS = (:production_hybrid, :memoized_dense, :independent_oracle)

@inline function _arg(values, name, default=nothing)
    index = findfirst(==(name), values)
    index === nothing && return default
    index == length(values) && throw(ArgumentError("missing value for $(name)"))
    return values[index + 1]
end

function _config(args)
    method = Symbol(_arg(args, "--method", "production_hybrid"))
    method in HYBRID_METHODS || throw(ArgumentError("unsupported hybrid shadow method $(method)"))
    xi = parse(Float64, _arg(args, "--xi", "0.0"))
    any(isapprox(xi, value; atol=1e-9, rtol=0.0) for value in HYBRID_XI) ||
        throw(ArgumentError("shadow xi must be one of -0.5, 0.0, 0.5"))
    calculation_sha = String(_arg(args, "--calculation-sha", ""))
    occursin(r"^[0-9a-fA-F]{40}$", calculation_sha) ||
        throw(ArgumentError("calculation-sha must be an immutable 40-character SHA"))
    output_dir = abspath(String(_arg(args, "--output-dir", joinpath(pwd(), "hybrid_shadow_artifact"))))
    tag = String(_arg(args, "--tag", "cep_cascade_production_shadow_v2"))
    return (; method, xi, calculation_sha, output_dir, tag)
end

@inline function _method_steps(method::Symbol)
    method === :production_hybrid && return (coarse=0.05, fine=0.025, policy=:rho_support_hybrid, levels=4)
    method === :memoized_dense && return (coarse=0.0125, fine=0.00625, policy=:uniform_nested, levels=1)
    return (coarse=0.00625, fine=0.003125, policy=:uniform_nested, levels=1)
end

function _rho_grid(step::Float64)
    count = Int(round(4.0 / step))
    return Float64.(collect(range(0.0, 4.0; length=count + 1)))
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

function _curve_rows(path::String, cfg, steps)
    isfile(path) || return NamedTuple[]
    rows = NamedTuple[]
    seen = Set{Tuple{Float64, Float64, Float64}}()
    for row in CSV.File(path)
        T = _safe_float(getproperty(row, :T_MeV))
        rho = _safe_float(getproperty(row, :rho))
        xi = _safe_float(getproperty(row, :xi))
        isfinite(T) && isfinite(rho) && isfinite(xi) || continue
        key = (T, xi, rho)
        key in seen && continue
        push!(seen, key)
        on_coarse = isapprox(rho / steps.coarse, round(rho / steps.coarse); atol=1e-8, rtol=0.0)
        on_fine = isapprox(rho / steps.fine, round(rho / steps.fine); atol=1e-8, rtol=0.0)
        role = on_coarse ? "coarse_grid" : on_fine ? "fine_grid" :
            cfg.method === :production_hybrid ? "stage_c_support" : "targeted"
        level = on_coarse ? 0 : on_fine ? 1 : cfg.method === :production_hybrid ? 4 : 2
        converged = _bool(getproperty(row, :converged))
        residual = _safe_float(getproperty(row, :residual_norm))
        push!(rows, (
            xi=xi,
            method=String(cfg.method),
            calculation_sha=cfg.calculation_sha,
            T_MeV=T,
            rho_level=level,
            rho=rho,
            muq_MeV=_safe_float(getproperty(row, :mu_avg_MeV)),
            pressure_fm4=_safe_float(getproperty(row, :pressure_fm4)),
            residual_norm=residual,
            iterations=Int(round(_safe_float(getproperty(row, :iterations)))),
            converged=converged,
            finite=converged && isfinite(residual),
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

function _anchor_record(diag, T::Float64)
    records = get(diag, "sweep_records", NamedTuple[])
    for record in records
        candidate = try Float64(getproperty(record, :T_MeV)) catch; NaN end
        isfinite(candidate) && isapprox(candidate, T; atol=1e-8, rtol=0.0) && return record
    end
    return nothing
end

@inline function _field(record, name::Symbol, default=nothing)
    record === nothing && return default
    try
        return getproperty(record, name)
    catch
        return default
    end
end

function _run_anchor(cfg, T::Float64, anchor_dir::String, steps, telemetry)
    mkpath(anchor_dir)
    result = Models.run_production_phase_pipeline(
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
        work_telemetry=telemetry,
        memoize_uniform=steps.policy === :uniform_nested,
        promote_reference=false,
    )
    return result
end

function _slice_row(cfg, T, result)
    diag = result.diagnostics
    record = _anchor_record(diag, T)
    status = String(_field(record, :status, result.cep.result_status))
    cache = get(diag, "rho_support_cache", nothing)
    cache_dict = cache isa AbstractDict ? cache : Dict{String, Any}()
    solver_failure_count = Int(_field(record, :stats_failure, get(diag, "scan_failure", 0)))
    return (
        xi=cfg.xi, method=String(cfg.method), calculation_sha=cfg.calculation_sha,
        T_MeV=T,
        stage_a_status=String(_field(record, :hybrid_stage_a_status, _field(record, :cascade_status, :not_run))),
        stage_b_status=String(_field(record, :hybrid_stage_b_status, :not_run)),
        stage_c_status=String(_field(record, :hybrid_stage_c_status, :not_run)),
        stage_used=String(_field(record, :stage_used, :uniform_nested)),
        upgrade_reason=String(_field(record, :hybrid_upgrade_reason, "not_applicable")),
        result_status=status,
        raw_status=String(_field(record, :raw_status, status)),
        geometry_converged=Bool(_field(record, :geometry_converged, false)),
        position_error_MeV=_safe_float(_field(record, :position_error_MeV, NaN)),
        density_error=_safe_float(_field(record, :density_error, NaN)),
        maxwell_area_gate=_safe_float(_field(record, :maxwell_area_gate, NaN)),
        area_residual=_safe_float(_field(record, :area_residual, NaN)),
        rho_hadron=_safe_float(_field(record, :rho_hadron, NaN)),
        rho_quark=_safe_float(_field(record, :rho_quark, NaN)),
        mu_transition_MeV=_safe_float(_field(record, :mu_transition, NaN)),
        mu_spinodal_hadron_MeV=_safe_float(_field(record, :mu_spinodal_hadron_MeV, NaN)),
        mu_spinodal_quark_MeV=_safe_float(_field(record, :mu_spinodal_quark_MeV, NaN)),
        rho_spinodal_hadron=_safe_float(_field(record, :rho_spinodal_hadron, NaN)),
        rho_spinodal_quark=_safe_float(_field(record, :rho_spinodal_quark, NaN)),
        support_low=_safe_float(_field(record, :hybrid_support_low, NaN)),
        support_high=_safe_float(_field(record, :hybrid_support_high, NaN)),
        support_point_count=Int(_field(record, :hybrid_verification_point_count, 0)),
        targeted_additions=Int(_field(record, :cascade_targeted_count, get(cache_dict, "targeted_additions", 0))),
        unique_solves=Int(_field(_field(record, :cache_stats, nothing), :unique_solves, get(cache_dict, "unique_solves", 0))),
        point_requests=Int(_field(_field(record, :cache_stats, nothing), :point_requests, get(cache_dict, "point_requests", 0))),
        cache_hits=Int(_field(_field(record, :cache_stats, nothing), :cache_hits, get(cache_dict, "cache_hits", 0))),
        finite_and_converged=solver_failure_count == 0,
        solver_failure_count=solver_failure_count,
    )
end

function _run_job(cfg)
    steps = _method_steps(cfg.method)
    started = time_ns()
    telemetry = Models.SolverWorkTelemetry()
    curve_rows = NamedTuple[]
    slice_rows = NamedTuple[]
    accuracy_rows = NamedTuple[]
    results = Models.PhasePipelineResult[]

    for T in HYBRID_ANCHORS[cfg.xi]
        anchor_dir = joinpath(cfg.output_dir, "anchors", "T_$(replace(string(T), "." => "p"))")
        result = _run_anchor(cfg, T, anchor_dir, steps, telemetry)
        push!(results, result)
        append!(curve_rows, _curve_rows(joinpath(anchor_dir, "trho_scan.csv"), cfg, steps))
        row = _slice_row(cfg, T, result)
        push!(slice_rows, row)
        cep = result.cep
        push!(accuracy_rows, (
            xi=cfg.xi, method=String(cfg.method), calculation_sha=cfg.calculation_sha,
            anchor_T_MeV=T, result_status=String(cep.result_status), found=cep.found,
            T_last_first_order_MeV=cep.T_last_first_order_MeV,
            muq_last_first_order_MeV=cep.mu_last_first_order_MeV,
            T_first_monotone_MeV=cep.T_first_monotone_MeV,
            ambiguity_width_T_MeV=cep.ambiguity_width_T_MeV,
            temperature_resolution_target_MeV=cep.temperature_resolution_target_MeV,
        ))
    end

    snapshot = Models.solver_work_snapshot(telemetry)
    elapsed = (time_ns() - started) / 1e9
    cache = Dict(
        "point_requests" => sum(Int(row.point_requests) for row in slice_rows),
        "cache_hits" => sum(Int(row.cache_hits) for row in slice_rows),
        "unique_solves" => sum(Int(row.unique_solves) for row in slice_rows),
        "targeted_additions" => sum(Int(row.targeted_additions) for row in slice_rows),
        "failed_points" => sum(Int(row.solver_failure_count) for row in slice_rows),
        "aggregation_scope" => "all_anchors",
    )
    cache["point_request_reconciliation"] = cache["point_requests"] == cache["unique_solves"] + cache["cache_hits"]
    cost = [(
        xi=cfg.xi, method=String(cfg.method), calculation_sha=cfg.calculation_sha,
        equilibrium_requests=snapshot.equilibrium_requests,
        fixedrho_requests=snapshot.fixedrho_requests,
        unique_solves=cache["unique_solves"], requested_point_calls=cache["point_requests"],
        uncached_equivalent_requests=cache["point_requests"], cache_hits=cache["cache_hits"],
        targeted_additions=cache["targeted_additions"],
        residual_calls=snapshot.nlsolve_f_calls + snapshot.postprocess_residual_calls,
        jacobian_calls=snapshot.nlsolve_g_calls,
        newton_iterations=snapshot.newton_iterations,
        trust_region_iterations=snapshot.trust_region_iterations,
        nonconverged_attempts=snapshot.nonconverged_attempts,
        fallback_count=snapshot.root_fallbacks,
        governed_rescue_count=snapshot.governed_rescues,
        retry_count=snapshot.scan_retries, exception_count=snapshot.exceptions,
        runner_seconds=elapsed,
    )]
    _write_csv(joinpath(cfg.output_dir, "curve_points.csv"), curve_rows)
    _write_csv(joinpath(cfg.output_dir, "slice_metrics.csv"), slice_rows)
    _write_csv(joinpath(cfg.output_dir, "method_costs.csv"), cost)
    _write_csv(joinpath(cfg.output_dir, "cep_accuracy.csv"), accuracy_rows)

    hashes = Dict{String, String}()
    for name in ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv")
        path = joinpath(cfg.output_dir, name)
        hashes[name] = bytes2hex(sha256(read(path)))
    end
    summary = Dict(
        "schema_version" => "cep_cascade_production_shadow_v2",
        "xi" => cfg.xi, "method" => String(cfg.method), "tag" => cfg.tag,
        "calculation_sha" => cfg.calculation_sha,
        "workflow_head_sha" => cfg.calculation_sha,
        "anchors" => HYBRID_ANCHORS[cfg.xi],
        "parameters" => Dict(
            "p_num" => 24, "t_num" => 8,
            "thermo_quadrature_policy" => "rs_reduced_adaptive",
            "thermo_quadrature_rtol" => 1e-8, "thermo_quadrature_atol" => 1e-10,
            "thermo_quadrature_maxevals" => 10^7, "iterations" => 80,
            "rho_refinement_policy" => String(steps.policy),
            "rho_coarse_step" => steps.coarse, "rho_fine_step" => steps.fine,
            "rho_refine_levels" => steps.levels,
            "rho_hybrid_local_step" => 0.003125,
            "rho_hybrid_support_padding" => 0.025,
            "rho_support_targeted_cap" => steps.policy === :rho_support_hybrid ? 12 : 0,
        ),
        "telemetry" => Dict(string(field) => getproperty(snapshot, field) for field in propertynames(snapshot)),
        "cache" => cache,
        "runner_seconds" => elapsed,
        "finite_and_converged_final" => all(row -> row.finite_and_converged, slice_rows),
        "curve_file_sha256" => hashes,
        "provenance" => Dict("calculation_sha" => cfg.calculation_sha,
            "reference_write" => false, "anchor_run_count" => length(HYBRID_ANCHORS[cfg.xi])),
    )
    open(joinpath(cfg.output_dir, "job_summary.json"), "w") do io
        JSON3.write(io, summary)
        write(io, '\n')
    end
    open(joinpath(cfg.output_dir, "manifest.json"), "w") do io
        JSON3.pretty(io, Dict("schema_version" => "cep_cascade_production_shadow_v2",
            "calculation_sha" => cfg.calculation_sha, "files" => hashes,
            "anchors" => HYBRID_ANCHORS[cfg.xi], "method" => String(cfg.method)))
        write(io, '\n')
    end
    println(JSON3.write(summary))
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    _run_job(_config(ARGS))
end
