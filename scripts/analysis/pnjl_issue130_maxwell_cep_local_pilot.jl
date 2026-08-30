#!/usr/bin/env julia

"""Run one fixed-(xi,T) Maxwell CEP-local diagnostic target.

The runner is intentionally separate from the production phase pipeline.  It
materializes the complete rho curve at the preflight step, then adds only
midpoints in the currently active three-crossing bracket.  Every refinement
level reruns the public strict candidate replay; no oracle labels or reference
artifacts are read or written.
"""

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using CSV
using JSON3
using SHA

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
end
const MODELS = Main.Models

include(joinpath(@__DIR__, "pnjl_maxwell_endpoint_candidate_feasibility.jl"))
const REPLAY = Main.MaxwellEndpointCandidateFeasibility

const DEFAULT_SCHEMA_VERSION = "pnjl_issue130_maxwell_cep_local_pilot_v1"
const BASE_STEP = 0.00625
const RHO_MAX = 4.0
const MAX_TARGETED = 12
const LOCAL_STEPS = (0.003125, 0.0015625, 0.00078125)
const POSITION_TOL = 0.025
const DENSITY_TOL = 0.0025
const AREA_GATE = 5e-5

function _arg(args, name, default=nothing)
    index = findfirst(==(name), args)
    index === nothing && return default
    index == length(args) && throw(ArgumentError("missing value for $name"))
    return args[index + 1]
end

function _config(args)
    calculation_sha = String(_arg(args, "--calculation-sha", ""))
    workflow_head_sha = String(_arg(args, "--workflow-head-sha", ""))
    target_id = String(_arg(args, "--target-id", ""))
    target_list = abspath(String(_arg(args, "--target-list", "")))
    output_dir = abspath(String(_arg(args, "--output-dir", joinpath(pwd(), "maxwell_cep_local_artifact"))))
    tag = String(_arg(args, "--tag", "issue130_maxwell_cep_local_pilot_v1"))
    selection = String(_arg(args, "--selection", "pilot_candidate"))
    schema_version = String(_arg(args, "--schema-version", DEFAULT_SCHEMA_VERSION))
    occursin(r"^[0-9a-fA-F]{40}$", calculation_sha) ||
        throw(ArgumentError("calculation-sha must be an immutable 40-character SHA"))
    occursin(r"^[0-9a-fA-F]{40}$", workflow_head_sha) ||
        throw(ArgumentError("workflow-head-sha must be an immutable 40-character SHA"))
    isempty(target_id) && throw(ArgumentError("target-id is required"))
    isfile(target_list) || throw(ArgumentError("target-list does not exist: $target_list"))
    isempty(selection) && throw(ArgumentError("selection must not be empty"))
    isempty(schema_version) && throw(ArgumentError("schema-version must not be empty"))
    return (; calculation_sha, workflow_head_sha, target_id, target_list, output_dir, tag, selection, schema_version)
end

function _field(row, name::Symbol, default=nothing)
    hasproperty(row, name) || return default
    value = getproperty(row, name)
    value === missing ? default : value
end

function _float(value, label)
    value === nothing && throw(ArgumentError("missing numeric target field: $label"))
    result = try Float64(value) catch; NaN end
    isfinite(result) || throw(ArgumentError("non-finite target field: $label"))
    result
end

function _target_row(path::String, target_id::String, selection::String)
    rows = collect(CSV.File(path))
    matches = filter(row -> String(_field(row, :target_id, "")) == target_id, rows)
    length(matches) == 1 || throw(ArgumentError("target-id must select exactly one row: $target_id"))
    row = first(matches)
    String(_field(row, :target_kind, "")) == "maxwell_fixed_xi_T" ||
        throw(ArgumentError("target is not a fixed-(xi,T) Maxwell target: $target_id"))
    String(_field(row, :pilot_selection, "")) == selection ||
        throw(ArgumentError("target does not match requested selection '$selection': $target_id"))
    (; target_id,
        xi=_float(_field(row, :xi), :xi),
        T_MeV=_float(_field(row, :T_MeV), :T_MeV),
        cep_T_low_MeV=_float(_field(row, :cep_T_low_MeV), :cep_T_low_MeV),
        cep_T_high_MeV=_float(_field(row, :cep_T_high_MeV), :cep_T_high_MeV),
        grid_status=String(_field(row, :grid_status, "")),
        reason=String(_field(row, :reason, "")))
end

function _rho_grid()
    count = Int(round(RHO_MAX / BASE_STEP))
    Float64.(collect(range(0.0, RHO_MAX; length=count + 1)))
end

function _curve(session, target)
    rows = collect(values(session.cache))
    sort!(rows; by=row -> row.rho)
    points = [REPLAY.CurvePoint(target.xi, target.T_MeV, "maxwell_cep_local_pilot",
        Float64(row.rho), Float64(row.muq_MeV), Float64(row.residual_norm),
        row.targeted ? 1 : 0, row.targeted ? "targeted_midpoint" : "base_grid")
        for row in rows if row.converged && row.finite]
    isempty(points) && return nothing
    REPLAY.CurveData(points, Float64[point.rho for point in points],
        Float64[point.mu for point in points])
end

function _bracket(curve, result)
    result.status == :first_order || return nothing
    crossings = result.candidate.crossings
    isempty(crossings) && return nothing
    crossing = first(crossings)
    index = clamp(searchsortedlast(curve.rho, crossing), 1, length(curve.rho) - 1)
    left, right = curve.rho[index], curve.rho[index + 1]
    (; left, right, width=right - left, positive=left > 0.0 && right > 0.0)
end

function _metric_row(target, level, targeted_count, session, curve, result, previous)
    candidate = result.status == :first_order ? result.candidate : nothing
    previous_candidate = previous === nothing ? nothing : previous.candidate
    position_error = candidate === nothing || previous_candidate === nothing ? NaN :
        abs(candidate.mu - previous_candidate.mu)
    density_error = if candidate === nothing || previous_candidate === nothing
        NaN
    else
        maximum(abs.(Float64[first(candidate.crossings), last(candidate.crossings)] .-
            Float64[first(previous_candidate.crossings), last(previous_candidate.crossings)]))
    end
    area_residual = candidate === nothing ? NaN : abs(Float64(candidate.area))
    geometry_evaluable = candidate !== nothing && previous_candidate !== nothing &&
        isfinite(position_error) && isfinite(density_error)
    geometry_converged = geometry_evaluable && position_error <= POSITION_TOL &&
        density_error <= DENSITY_TOL && area_residual <= AREA_GATE
    bracket = curve === nothing ? nothing : _bracket(curve, result)
    snapshot = MODELS.TrhoScan.rho_session_snapshot(session)
    (
        target_id=target.target_id, xi=target.xi, T_MeV=target.T_MeV,
        level=level, targeted_additions=targeted_count,
        status=String(result.status), reason=String(result.reason),
        candidate_count=length(result.roots),
        candidate_mu=candidate === nothing ? missing : candidate.mu,
        candidate_area=candidate === nothing ? missing : candidate.area,
        rho_hadron=candidate === nothing ? missing : first(candidate.crossings),
        rho_quark=candidate === nothing ? missing : last(candidate.crossings),
        endpoint_dependent=get(result, :endpoint_dependent, false),
        bracket_low=bracket === nothing ? missing : bracket.left,
        bracket_high=bracket === nothing ? missing : bracket.right,
        bracket_width=bracket === nothing ? missing : bracket.width,
        positive_rho_bracket=bracket === nothing ? false : bracket.positive,
        position_error_MeV=position_error, density_error=density_error,
        area_residual=area_residual, geometry_evaluable=geometry_evaluable,
        geometry_converged=geometry_converged,
        point_requests=snapshot.point_requests, unique_solves=snapshot.unique_solves,
        cache_hits=snapshot.cache_hits, solver_failures=snapshot.failed_points,
    )
end

function _target_for_next(curve, result, selected::Set{Float64})
    bracket = _bracket(curve, result)
    bracket === nothing && return nothing
    target = 0.5 * (bracket.left + bracket.right)
    target in selected && return nothing
    target
end

function _write_csv(path, rows)
    mkpath(dirname(path))
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _write_json(path, value)
    open(path, "w") do io
        JSON3.pretty(io, value)
        write(io, '\n')
    end
end

function _write_job(cfg)
    target = _target_row(cfg.target_list, cfg.target_id, cfg.selection)
    mkpath(cfg.output_dir)
    telemetry = MODELS.SolverWorkTelemetry()
    session = MODELS.TrhoScan.new_rho_point_session(
        model_kind=:PNJL, reverse_rho=true, seed_policy=:candidates,
        solver_backend=:models, p_num=24, t_num=8, iterations=80,
        thermo_quadrature_policy=:rs_reduced_adaptive,
        thermo_quadrature_rtol=1e-8, thermo_quadrature_atol=1e-10,
        thermo_quadrature_maxevals=10^7, telemetry=telemetry)

    for rho in reverse(_rho_grid())
        MODELS.TrhoScan.rho_point!(session, target.T_MeV, target.xi, rho)
    end

    curve = _curve(session, target)
    result = curve === nothing ? (status=:invalid, reason="solver_or_curve_failure",
        roots=NamedTuple[], candidate=nothing) : REPLAY.strict_candidate(curve)
    metrics = NamedTuple[_metric_row(target, 0, 0, session, curve, result, nothing)]
    previous = result
    selected = Set{Float64}()
    while length(selected) < MAX_TARGETED
        curve === nothing && break
        target_rho = _target_for_next(curve, result, selected)
        target_rho === nothing && break
        push!(selected, target_rho)
        MODELS.TrhoScan.rho_point!(session, target.T_MeV, target.xi, target_rho; targeted=true)
        curve = _curve(session, target)
        result = curve === nothing ? (status=:invalid, reason="solver_or_curve_failure",
            roots=NamedTuple[], candidate=nothing) : REPLAY.strict_candidate(curve)
        push!(metrics, _metric_row(target, length(selected), length(selected), session,
            curve, result, previous))
        previous = result
    end

    final_metric = last(metrics)
    final_candidate = final_metric.candidate_count == 1 && final_metric.status == "first_order"
    geometry_pass = any(row -> row.geometry_converged, metrics)
    snapshot = MODELS.solver_work_snapshot(telemetry)
    curve_rows = NamedTuple[]
    if curve !== nothing
        for point in curve.points
            push!(curve_rows, (
                target_id=target.target_id, xi=point.xi, T_MeV=point.T,
                rho=point.rho, muq_MeV=point.mu, residual_norm=point.residual,
                rho_level=point.rho_level, sampling_role=point.sampling_role,
                converged=true, finite=true,
                calculation_sha=cfg.calculation_sha,
                workflow_head_sha=cfg.workflow_head_sha,
            ))
        end
    end
    _write_csv(joinpath(cfg.output_dir, "curve_points.csv"), curve_rows)
    _write_csv(joinpath(cfg.output_dir, "slice_metrics.csv"), metrics)
    _write_csv(joinpath(cfg.output_dir, "policy_frontier.csv"), [(
        target_id=target.target_id, cap=MAX_TARGETED,
        selected_points=length(selected), candidate_unique=final_candidate,
        geometry_converged=geometry_pass,
        policy="cross_level_candidate_stability_v1",
    )])
    _write_csv(joinpath(cfg.output_dir, "method_costs.csv"), [(
        target_id=target.target_id, xi=target.xi, T_MeV=target.T_MeV,
        calculation_sha=cfg.calculation_sha, workflow_head_sha=cfg.workflow_head_sha,
        base_rho_step=BASE_STEP, base_rho_max=RHO_MAX,
        equilibrium_requests=snapshot.equilibrium_requests,
        fixedrho_requests=snapshot.fixedrho_requests,
        unique_solves=snapshot.fixedrho_requests,
        point_requests=session.point_requests, cache_hits=session.cache_hits,
        targeted_additions=session.targeted_additions,
        residual_calls=snapshot.nlsolve_f_calls + snapshot.postprocess_residual_calls,
        jacobian_calls=snapshot.nlsolve_g_calls,
        newton_iterations=snapshot.newton_iterations,
        nonconverged_attempts=snapshot.nonconverged_attempts,
        fallback_count=snapshot.root_fallbacks, retry_count=snapshot.scan_retries,
        exception_count=snapshot.exceptions, solver_failures=session.failed_points,
    )])

    verdict = if session.failed_points > 0 || curve === nothing
        "solver_or_curve_failure"
    elseif final_candidate && geometry_pass
        "candidate_and_geometry_feasible"
    elseif result.status == :multiple_candidates
        "multiple_maxwell_candidate_inconclusive"
    else
        "candidate_or_geometry_inconclusive"
    end
    summary = Dict(
        "schema_version" => cfg.schema_version, "target_id" => target.target_id,
        "selection" => cfg.selection,
        "xi" => target.xi, "T_MeV" => target.T_MeV,
        "grid_status" => target.grid_status, "preflight_reason" => target.reason,
        "verdict" => verdict, "final_status" => final_metric.status,
        "final_reason" => final_metric.reason,
        "final_candidate_count" => final_metric.candidate_count,
        "final_candidate_mu_MeV" => final_metric.candidate_mu,
        "final_area_residual" => final_metric.area_residual,
        "final_geometry_converged" => geometry_pass,
        "targeted_additions" => length(selected),
        "curve_points" => length(curve_rows),
        "finite_and_converged" => curve !== nothing && session.failed_points == 0,
        "reference_write" => false, "oracle_labels_consumed" => false,
        "solver_called" => true,
    )
    _write_json(joinpath(cfg.output_dir, "target_summary.json"), summary)
    _write_json(joinpath(cfg.output_dir, "provenance.json"), Dict(
        "schema_version" => cfg.schema_version, "target_id" => target.target_id,
        "selection" => cfg.selection,
        "calculation_sha" => cfg.calculation_sha,
        "workflow_head_sha" => cfg.workflow_head_sha,
        "target_list" => cfg.target_list, "tag" => cfg.tag,
        "solver_called" => true, "reference_write" => false,
        "oracle_labels_consumed" => false,
        "base_rho_step" => BASE_STEP, "base_rho_max" => RHO_MAX,
        "targeted_cap" => MAX_TARGETED,
    ))
    files = Dict{String,String}()
    for (root, _, names) in walkdir(cfg.output_dir)
        for name in names
            name == "manifest.json" && continue
            path = joinpath(root, name)
            files[replace(relpath(path, cfg.output_dir), '\\' => '/')] = bytes2hex(open(sha256, path))
        end
    end
    manifest = Dict(
        "schema_version" => cfg.schema_version, "verdict" => verdict,
        "target_id" => target.target_id, "xi" => target.xi,
        "T_MeV" => target.T_MeV, "calculation_sha" => cfg.calculation_sha,
        "workflow_head_sha" => cfg.workflow_head_sha, "tag" => cfg.tag,
        "selection" => cfg.selection,
        "reference_write" => false, "oracle_labels_consumed" => false,
        "solver_called" => true, "targeted_cap" => MAX_TARGETED,
        "files" => files, "telemetry" => Dict(string(field) => getproperty(snapshot, field)
            for field in propertynames(snapshot)),
    )
    _write_json(joinpath(cfg.output_dir, "manifest.json"), manifest)
    println(JSON3.write(manifest))
end

if abspath(PROGRAM_FILE) == @__FILE__
    _write_job(_config(ARGS))
end
