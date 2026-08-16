#!/usr/bin/env julia

"""Run one targeted Issue #130 closure job.

The production shadow runner is included so the numerical path stays exactly
the calculation-SHA path.  This wrapper only supplies the frozen target list
and writes a small, independently versioned artifact containing raw rho-mu
curves, slice diagnostics, and solver-work telemetry.  It never writes a
reference artifact.
"""

include(joinpath(@__DIR__, "pnjl_cep_hybrid_production_shadow.jl"))

const TARGETED_CLOSURE_SCHEMA = "pnjl_c2_targeted_closure_job_v1"
const TARGETED_CLOSURE_CALCULATION_SHA =
    "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"

const REGRESSION_TARGETS = [
    (id="xi_m0p35_T51", xi=-0.35, T_MeV=51.0),
    (id="xi_m0p25_T41", xi=-0.25, T_MeV=41.0),
    (id="xi_m0p20_T41", xi=-0.20, T_MeV=41.0),
    (id="xi_m0p15_T41", xi=-0.15, T_MeV=41.0),
    (id="xi_m0p10_T41", xi=-0.10, T_MeV=41.0),
    (id="xi_0p00_T51", xi=0.0, T_MeV=51.0),
    (id="xi_0p30_T21", xi=0.30, T_MeV=21.0),
    (id="xi_0p35_T51", xi=0.35, T_MeV=51.0),
    (id="xi_0p35_T101", xi=0.35, T_MeV=101.0),
]

const CEP_TARGETS = [
    (id="xi_0p125_CEP_midpoint", xi=0.125, T_MeV=126.25,
        bracket_low=126.1875, bracket_high=126.3125),
    (id="xi_0p39375_CEP_midpoint", xi=0.39375, T_MeV=113.5,
        bracket_low=113.4375, bracket_high=113.5625),
    (id="xi_0p50_CEP_midpoint", xi=0.5, T_MeV=107.0,
        bracket_low=106.9375, bracket_high=107.0625),
]

@inline function _target_arg(args, name, default=nothing)
    index = findfirst(==(name), args)
    index === nothing && return default
    index == length(args) && throw(ArgumentError("missing value for $(name)"))
    args[index + 1]
end

function _target_config(args)
    scope = Symbol(_target_arg(args, "--scope", "regression_curves"))
    scope in (:regression_curves, :cep_brackets) ||
        throw(ArgumentError("numerical targeted jobs support regression_curves or cep_brackets"))
    target_id = String(_target_arg(args, "--target-id", ""))
    targets = scope === :regression_curves ? REGRESSION_TARGETS : CEP_TARGETS
    target = findfirst(item -> item.id == target_id, targets)
    target === nothing && throw(ArgumentError("unknown target-id $(target_id) for $(scope)"))
    calculation_sha = lowercase(String(_target_arg(args, "--calculation-sha", "")))
    occursin(r"^[0-9a-f]{40}$", calculation_sha) ||
        throw(ArgumentError("calculation-sha must be a lowercase 40-character SHA"))
    calculation_sha == TARGETED_CLOSURE_CALCULATION_SHA ||
        throw(ArgumentError("calculation SHA does not match frozen targeted-closure input"))
    methods = Tuple(Symbol.(split(String(_target_arg(args, "--methods",
        "production_hybrid,independent_oracle")), ",")))
    all(method -> method in (:production_hybrid, :independent_oracle), methods) ||
        throw(ArgumentError("methods must be production_hybrid and/or independent_oracle"))
    isempty(methods) && throw(ArgumentError("at least one method is required"))
    output_dir = abspath(String(_target_arg(args, "--output-dir", joinpath(pwd(), "targeted_job"))))
    tag = String(_target_arg(args, "--tag", "issue130_c2_targeted_closure_v1_20260816"))
    postprocess_sha = lowercase(String(_target_arg(args, "--postprocess-sha", "unknown")))
    source_run_id = String(_target_arg(args, "--source-run-id", "pending"))
    (; scope, target_id, target, methods, calculation_sha, output_dir, tag,
        postprocess_sha, source_run_id)
end

function _target_properties(target)
    properties = Dict{String, Any}(
        "id" => target.id, "xi" => target.xi, "T_MeV" => target.T_MeV)
    hasproperty(target, :bracket_low) && (properties["bracket_low_MeV"] = target.bracket_low)
    hasproperty(target, :bracket_high) && (properties["bracket_high_MeV"] = target.bracket_high)
    properties
end

function _hash_file(path)
    bytes2hex(SHA.sha256(read(path)))
end

function _run_method(config)
    method_root = joinpath(config.output_dir, String(config.method))
    mkpath(method_root)
    steps = _method_steps(config.method)
    telemetry = Models.SolverWorkTelemetry()
    started = time_ns()
    curve_rows = NamedTuple[]
    slice_rows = NamedTuple[]
    accuracy_rows = NamedTuple[]

    target = config.target
    anchor_dir = joinpath(method_root, "anchor")
    result = _run_anchor((;
        method=config.method, xi=target.xi, calculation_sha=config.calculation_sha,
        tag=config.tag, scope=config.scope,
        endpoint_policy=:three_crossing_endpoint_local_v2,
        candidate_policy=:unique_three_crossing_sign_change_v2),
        target.T_MeV, anchor_dir, steps, telemetry)
    anchor_cfg = (;
        method=config.method, xi=target.xi, calculation_sha=config.calculation_sha,
        tag=config.tag, scope=config.scope,
        endpoint_policy=:three_crossing_endpoint_local_v2,
        candidate_policy=:unique_three_crossing_sign_change_v2)
    curves = _curve_rows(joinpath(anchor_dir, "trho_scan.csv"), anchor_cfg, steps)
    append!(curve_rows, curves)
    slice = _slice_row(anchor_cfg, target.T_MeV, result)
    push!(slice_rows, slice)
    push!(accuracy_rows, (
        target_id=target.id, xi=target.xi, method=String(config.method),
        T_MeV=target.T_MeV, result_status=String(result.cep.result_status),
        found=result.cep.found,
        T_last_first_order_MeV=result.cep.T_last_first_order_MeV,
        T_first_monotone_MeV=result.cep.T_first_monotone_MeV,
        ambiguity_width_T_MeV=result.cep.ambiguity_width_T_MeV,
        hybrid_status=slice.result_status,
        maxwell_candidate_count=slice.maxwell_candidate_count,
        maxwell_crossing_count=slice.maxwell_crossing_count,
        geometry_converged=slice.geometry_converged,
        position_error_MeV=slice.position_error_MeV,
        density_error=slice.density_error,
        maxwell_area_gate=slice.maxwell_area_gate,
        area_residual=slice.area_residual,
        finite_and_converged=slice.finite_and_converged,
    ))

    snapshot = Models.solver_work_snapshot(telemetry)
    elapsed = (time_ns() - started) / 1e9
    cache = Dict(
        "point_requests" => Int(slice.point_requests),
        "cache_hits" => Int(slice.cache_hits),
        "unique_solves" => Int(slice.unique_solves),
        "targeted_additions" => Int(slice.targeted_additions),
        "failed_points" => Int(slice.solver_failure_count),
    )
    cache["point_request_reconciliation"] =
        cache["point_requests"] == cache["unique_solves"] + cache["cache_hits"]
    costs = [(
        target_id=target.id, xi=target.xi, T_MeV=target.T_MeV,
        method=String(config.method), calculation_sha=config.calculation_sha,
        unique_solves=cache["unique_solves"], point_requests=cache["point_requests"],
        cache_hits=cache["cache_hits"], targeted_additions=cache["targeted_additions"],
        failed_points=cache["failed_points"],
        equilibrium_requests=snapshot.equilibrium_requests,
        fixedrho_requests=snapshot.fixedrho_requests,
        residual_calls=snapshot.nlsolve_f_calls + snapshot.postprocess_residual_calls,
        jacobian_calls=snapshot.nlsolve_g_calls,
        newton_iterations=snapshot.newton_iterations,
        fallback_count=snapshot.root_fallbacks, retry_count=snapshot.scan_retries,
        exception_count=snapshot.exceptions, runner_seconds=elapsed,
    )]

    _write_csv(joinpath(method_root, "curve_points.csv"), curve_rows)
    _write_csv(joinpath(method_root, "slice_diagnostics.csv"), slice_rows)
    _write_csv(joinpath(method_root, "accuracy.csv"), accuracy_rows)
    _write_csv(joinpath(method_root, "method_costs.csv"), costs)
    files = Dict(name => _hash_file(joinpath(method_root, name)) for name in
        ("curve_points.csv", "slice_diagnostics.csv", "accuracy.csv", "method_costs.csv"))
    summary = Dict(
        "schema_version" => TARGETED_CLOSURE_SCHEMA,
        "scope" => String(config.scope), "target_id" => target.id,
        "target" => _target_properties(target), "method" => String(config.method),
        "tag" => config.tag, "source_run_id" => config.source_run_id,
        "calculation_sha" => config.calculation_sha,
        "postprocess_sha" => config.postprocess_sha,
        "workflow_head_sha" => config.postprocess_sha,
        "rho_coarse_step" => steps.coarse, "rho_fine_step" => steps.fine,
        "rho_refinement_policy" => String(steps.policy),
        "endpoint_policy" => "three_crossing_endpoint_local_v2",
        "candidate_policy" => "unique_three_crossing_sign_change_v2",
        "solver_called" => true,
        "finite_and_converged_final" => all(row -> row.finite_and_converged, slice_rows),
        "cache" => cache,
        "telemetry" => Dict(string(field) => getproperty(snapshot, field)
            for field in propertynames(snapshot)),
        "runner_seconds" => elapsed, "files" => files,
        "provenance" => Dict(
            "reference_write" => false,
            "oracle_labels_used_for_routing" => false,
            "route_before_oracle_gate" => true,
            "targeted_overlay_only" => true,
        ),
    )
    open(joinpath(method_root, "manifest.json"), "w") do io
        JSON3.pretty(io, summary); write(io, '\n')
    end
    summary
end

function _run_target(config)
    mkpath(config.output_dir)
    summaries = Any[]
    for method in config.methods
        push!(summaries, _run_method((; config..., method=method)))
    end
    summary = Dict(
        "schema_version" => TARGETED_CLOSURE_SCHEMA,
        "scope" => String(config.scope), "target_id" => config.target_id,
        "target" => _target_properties(config.target), "tag" => config.tag,
        "source_run_id" => config.source_run_id,
        "calculation_sha" => config.calculation_sha,
        "postprocess_sha" => config.postprocess_sha,
        "methods" => String.(config.methods), "solver_called" => true,
        "method_summaries" => summaries,
        "provenance" => Dict(
            "reference_write" => false,
            "oracle_labels_used_for_routing" => false,
            "route_before_oracle_gate" => true,
            "targeted_overlay_only" => true,
        ),
    )
    open(joinpath(config.output_dir, "job_summary.json"), "w") do io
        JSON3.pretty(io, summary); write(io, '\n')
    end
    println(JSON3.write(summary))
end

if abspath(PROGRAM_FILE) == @__FILE__
    _run_target(_target_config(ARGS))
end
