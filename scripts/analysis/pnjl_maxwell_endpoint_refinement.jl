#!/usr/bin/env julia

"""Actions-only diagnostic for the low-rho Maxwell endpoint.

This runner solves one PNJL anchor, using the existing request-scoped rho
session.  It evaluates the complete 0:0.003125:4 curve first, then adds only
midpoints in the left outer-crossing bracket.  It is diagnostic-only and never
writes phase-reference artifacts.
"""

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using CSV
using JSON3
using SHA
using Statistics

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
end
const MODELS = Main.Models

include(joinpath(@__DIR__, "pnjl_maxwell_endpoint_candidate_feasibility.jl"))
const REPLAY = Main.MaxwellEndpointCandidateFeasibility

const XI = -0.5
const TEMPERATURE = 5.0
const RHO_MAX = 4.0
const BASE_STEP = 0.003125
const LOCAL_STEPS = (0.0015625, 0.00078125, 0.000390625, 0.0001953125)
const TARGETED_CAPS = (4, 6, 8, 10, 12)
const MAX_TARGETED = maximum(TARGETED_CAPS)
const POSITION_TOL = 0.025
const DENSITY_TOL = 0.0025
const AREA_TOL = 5e-6
const OUTER_AREA_GATE = 5e-5

function _arg(args, name, default=nothing)
    index = findfirst(==(name), args)
    index === nothing && return default
    index == length(args) && throw(ArgumentError("missing value for $name"))
    return args[index + 1]
end

function _config(args)
    calculation_sha = String(_arg(args, "--calculation-sha", ""))
    workflow_head_sha = String(_arg(args, "--workflow-head-sha", ""))
    occursin(r"^[0-9a-fA-F]{40}$", calculation_sha) ||
        throw(ArgumentError("calculation-sha must be an immutable 40-character SHA"))
    occursin(r"^[0-9a-fA-F]{40}$", workflow_head_sha) ||
        throw(ArgumentError("workflow-head-sha must be an immutable 40-character SHA"))
    output_dir = abspath(String(_arg(args, "--output-dir", joinpath(pwd(), "endpoint_refinement_artifact"))))
    tag = String(_arg(args, "--tag", "pnjl_maxwell_endpoint_refinement_v1"))
    return (; calculation_sha, workflow_head_sha, output_dir, tag)
end

function _rho_grid(step::Float64)
    count = Int(round(RHO_MAX / step))
    Float64.(collect(range(0.0, RHO_MAX; length=count + 1)))
end

function _point(curve_row, cfg, role::String, level::Int)
    REPLAY.CurvePoint(XI, TEMPERATURE, "endpoint_refinement", Float64(curve_row.rho),
        Float64(curve_row.muq_MeV), Float64(curve_row.residual_norm), level, role)
end

function _curve(session)
    rows = collect(values(session.cache))
    sort!(rows; by=row -> row.rho)
    points = [_point(row, nothing, row.targeted ? "targeted" : "base_grid",
        row.targeted ? 1 : 0) for row in rows if row.converged && row.finite]
    return REPLAY.CurveData(points, Float64[point.rho for point in points],
        Float64[point.mu for point in points])
end

function _bracket(curve, result)
    result.status == :first_order || return nothing
    crossings = result.candidate.crossings
    isempty(crossings) && return nothing
    crossing = first(crossings)
    index = searchsortedlast(curve.rho, crossing)
    index = clamp(index, 1, length(curve.rho) - 1)
    left, right = curve.rho[index], curve.rho[index + 1]
    return (left=left, right=right, width=right - left,
        positive=left > 0.0 && right > 0.0)
end

function _metric_row(level, targeted_count, session, curve, result, previous)
    snapshot = MODELS.TrhoScan.rho_session_snapshot(session)
    bracket = _bracket(curve, result)
    candidate = result.status == :first_order ? result.candidate : nothing
    position_error = if previous === nothing || candidate === nothing || previous.candidate === nothing
        NaN
    else
        abs(candidate.mu - previous.candidate.mu)
    end
    density_error = if previous === nothing || candidate === nothing || previous.candidate === nothing
        NaN
    else
        maximum(abs.(Float64[first(candidate.crossings), last(candidate.crossings)] .-
            Float64[first(previous.candidate.crossings), last(previous.candidate.crossings)]))
    end
    geometry = candidate !== nothing && previous !== nothing && previous.candidate !== nothing &&
        isfinite(position_error) && isfinite(density_error) &&
        position_error <= POSITION_TOL && density_error <= DENSITY_TOL &&
        abs(candidate.area) <= OUTER_AREA_GATE
    (
        level=level,
        targeted_additions=targeted_count,
        status=String(result.status),
        reason=String(result.reason),
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
        position_error_MeV=position_error,
        density_error=density_error,
        geometry_converged=geometry,
        point_requests=snapshot.point_requests,
        unique_solves=snapshot.unique_solves,
        cache_hits=snapshot.cache_hits,
        solver_failures=snapshot.failed_points,
    )
end

function _frontier(metrics)
    rows = NamedTuple[]
    for cap in TARGETED_CAPS
        subset = [row for row in metrics if row.targeted_additions <= cap]
        isempty(subset) && continue
        row = last(subset)
        selected_step = missing
        for step in LOCAL_STEPS
            if row.positive_rho_bracket && row.bracket_width <= step && row.geometry_converged
                selected_step = step
                break
            end
        end
        push!(rows, merge(row, (
            policy_cap=cap,
            selected_local_step=selected_step,
            policy_pass=selected_step !== missing && row.solver_failures == 0,
        )))
    end
    return rows
end

function _target_for_next(curve, result, target_set::Set{Float64})
    bracket = _bracket(curve, result)
    bracket === nothing && return nothing
    if bracket.width <= last(LOCAL_STEPS) && bracket.positive
        return nothing
    end
    target = 0.5 * (bracket.left + bracket.right)
    target in target_set && return nothing
    return target
end

function _write_csv(path, rows)
    mkpath(dirname(path))
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _write_job(cfg)
    mkpath(cfg.output_dir)
    telemetry = MODELS.SolverWorkTelemetry()
    session = MODELS.TrhoScan.new_rho_point_session(
        model_kind=:PNJL, reverse_rho=true, seed_policy=:candidates,
        solver_backend=:models, p_num=24, t_num=8, iterations=80,
        thermo_quadrature_policy=:rs_reduced_adaptive,
        thermo_quadrature_rtol=1e-8, thermo_quadrature_atol=1e-10,
        thermo_quadrature_maxevals=10^7, telemetry=telemetry)

    base = _rho_grid(BASE_STEP)
    for rho in reverse(base)
        MODELS.TrhoScan.rho_point!(session, TEMPERATURE, XI, rho)
    end

    metrics = NamedTuple[]
    targeted = Set{Float64}()
    curve = _curve(session)
    result = REPLAY.strict_candidate(curve)
    previous = nothing
    push!(metrics, _metric_row(0, 0, session, curve, result, previous))

    while length(targeted) < MAX_TARGETED
        target = _target_for_next(curve, result, targeted)
        target === nothing && break
        push!(targeted, target)
        MODELS.TrhoScan.rho_point!(session, TEMPERATURE, XI, target; targeted=true)
        curve = _curve(session)
        previous = result
        result = REPLAY.strict_candidate(curve)
        push!(metrics, _metric_row(length(targeted), length(targeted), session, curve, result, previous))
    end

    frontier = _frontier(metrics)
    selected = filter(row -> row.policy_pass, frontier)
    selected_policy = isempty(selected) ? nothing : first(sort(selected;
        by=row -> (row.selected_local_step, row.policy_cap)))
    snapshot = MODELS.solver_work_snapshot(telemetry)
    rows = NamedTuple[]
    for row in sort!(collect(values(session.cache)); by=row -> row.rho)
        push!(rows, (
            xi=XI, T_MeV=TEMPERATURE, rho=row.rho,
            muq_MeV=row.muq_MeV, pressure_fm4=row.pressure_fm4,
            residual_norm=row.residual_norm, iterations=row.result === nothing ? 0 : row.result.iterations,
            converged=row.converged, finite=row.finite,
            sampling_role=row.targeted ? "targeted_endpoint_local" : "base_grid",
            rho_level=row.targeted ? 1 : 0,
            calculation_sha=cfg.calculation_sha,
            workflow_head_sha=cfg.workflow_head_sha,
        ))
    end
    _write_csv(joinpath(cfg.output_dir, "curve_points.csv"), rows)
    _write_csv(joinpath(cfg.output_dir, "slice_metrics.csv"), metrics)
    _write_csv(joinpath(cfg.output_dir, "policy_frontier.csv"), frontier)
    _write_csv(joinpath(cfg.output_dir, "method_costs.csv"), [(
        xi=XI, T_MeV=TEMPERATURE, calculation_sha=cfg.calculation_sha,
        workflow_head_sha=cfg.workflow_head_sha,
        equilibrium_requests=snapshot.equilibrium_requests,
        fixedrho_requests=snapshot.fixedrho_requests,
        unique_solves=snapshot.fixedrho_requests,
        point_requests=session.point_requests, cache_hits=session.cache_hits,
        targeted_additions=session.targeted_additions,
        residual_calls=snapshot.nlsolve_f_calls + snapshot.postprocess_residual_calls,
        jacobian_calls=snapshot.nlsolve_g_calls,
        newton_iterations=snapshot.newton_iterations,
        trust_region_iterations=snapshot.trust_region_iterations,
        nonconverged_attempts=snapshot.nonconverged_attempts,
        fallback_count=snapshot.root_fallbacks,
        governed_rescue_count=snapshot.governed_rescues,
        retry_count=snapshot.scan_retries, exception_count=snapshot.exceptions,
    )])
    verdict = if selected_policy !== nothing
        "candidate_and_endpoint_feasible"
    elseif result.status == :multiple_candidates
        "multiple_maxwell_candidate_inconclusive"
    elseif session.failed_points > 0 || result.status == :invalid
        "solver_or_curve_failure"
    else
        "candidate_only_endpoint_inconclusive"
    end
    selected_value = selected_policy === nothing ? nothing : Dict(
        "cap" => selected_policy.policy_cap,
        "local_step" => selected_policy.selected_local_step,
        "targeted_additions" => selected_policy.targeted_additions,
    )
    files = Dict{String, String}()
    for (root, _, names) in walkdir(cfg.output_dir)
        for name in names
            name == "manifest.json" && continue
            path = joinpath(root, name)
            files[replace(relpath(path, cfg.output_dir), '\\' => '/')] = bytes2hex(open(sha256, path))
        end
    end
    _write_json = function(path, value)
        open(path, "w") do io
            JSON3.pretty(io, value)
            write(io, '\n')
        end
    end
    _write_json(joinpath(cfg.output_dir, "selected_policy.json"), Dict(
        "schema_version" => "pnjl_maxwell_endpoint_refinement_v1",
        "verdict" => verdict,
        "policy" => selected_value,
        "candidate_policy" => "unique_three_crossing_topology_v1",
        "endpoint_policy" => "positive_rho_bracket_required",
        "area_solver_tol" => AREA_TOL,
        "outer_area_gate" => OUTER_AREA_GATE,
        "position_tol_MeV" => POSITION_TOL,
        "density_tol" => DENSITY_TOL,
    ))
    manifest = Dict(
        "schema_version" => "pnjl_maxwell_endpoint_refinement_v1",
        "verdict" => verdict,
        "xi" => XI, "temperature_MeV" => TEMPERATURE,
        "base_rho_step" => BASE_STEP, "base_rho_max" => RHO_MAX,
        "local_step_candidates" => collect(LOCAL_STEPS),
        "targeted_cap_candidates" => collect(TARGETED_CAPS),
        "calculation_sha" => cfg.calculation_sha,
        "workflow_head_sha" => cfg.workflow_head_sha,
        "tag" => cfg.tag,
        "reference_write" => false,
        "files" => files,
        "telemetry" => Dict(string(field) => getproperty(snapshot, field) for field in propertynames(snapshot)),
        "selected_policy" => selected_value,
    )
    _write_json(joinpath(cfg.output_dir, "manifest.json"), manifest)
    println(JSON3.write(manifest))
end

if abspath(PROGRAM_FILE) == @__FILE__
    _write_job(_config(ARGS))
end
