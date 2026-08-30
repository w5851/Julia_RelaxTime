"""Request-scoped solver work counters.

The telemetry object is deliberately separate from `SolverResult v1`: callers
opt in by passing `work_telemetry=SolverWorkTelemetry()` to a solver request,
and can inspect a snapshot after the request or scan completes.  Keeping the
counters request-scoped avoids hidden global state and makes parallel pilot
jobs independently auditable.
"""

mutable struct SolverWorkTelemetry
    equilibrium_requests::Int
    fixedrho_requests::Int
    governed_attempts::Int
    nlsolve_f_calls::Int
    nlsolve_g_calls::Int
    postprocess_residual_calls::Int
    newton_iterations::Int
    trust_region_iterations::Int
    nonconverged_attempts::Int
    root_fallbacks::Int
    governed_rescues::Int
    exceptions::Int
    scan_retries::Int
end

SolverWorkTelemetry() = SolverWorkTelemetry(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
)

@inline function _require_solver_work_telemetry(telemetry)
    telemetry === nothing && return nothing
    telemetry isa SolverWorkTelemetry || throw(ArgumentError(
        "work_telemetry must be SolverWorkTelemetry or nothing, got $(typeof(telemetry))",
    ))
    return telemetry
end

@inline function record_solver_request!(telemetry; fixedrho::Bool=false)
    telemetry = _require_solver_work_telemetry(telemetry)
    telemetry === nothing && return nothing
    telemetry.equilibrium_requests += 1
    fixedrho && (telemetry.fixedrho_requests += 1)
    return nothing
end

@inline function record_governed_attempt!(telemetry; origin::Symbol=:primary)
    telemetry = _require_solver_work_telemetry(telemetry)
    telemetry === nothing && return nothing
    telemetry.governed_attempts += 1
    if origin === :method_rescue || origin === :seed_rescue || origin === :fallback
        telemetry.governed_rescues += 1
        origin === :method_rescue && (telemetry.root_fallbacks += 1)
    end
    return nothing
end

@inline function record_nlsolve_work!(telemetry, result, method::Symbol)
    telemetry = _require_solver_work_telemetry(telemetry)
    telemetry === nothing && return nothing
    telemetry.nlsolve_f_calls += Int(getproperty(result, :f_calls))
    telemetry.nlsolve_g_calls += Int(getproperty(result, :g_calls))
    iterations = Int(getproperty(result, :iterations))
    if method === :newton
        telemetry.newton_iterations += iterations
    elseif method === :trust_region
        telemetry.trust_region_iterations += iterations
    end
    return nothing
end

@inline function record_postprocess_residual!(telemetry)
    telemetry = _require_solver_work_telemetry(telemetry)
    telemetry === nothing && return nothing
    telemetry.postprocess_residual_calls += 1
    return nothing
end

@inline function record_attempt_outcome!(telemetry; converged::Bool)
    telemetry = _require_solver_work_telemetry(telemetry)
    telemetry === nothing && return nothing
    !converged && (telemetry.nonconverged_attempts += 1)
    return nothing
end

@inline function record_solver_exception!(telemetry)
    telemetry = _require_solver_work_telemetry(telemetry)
    telemetry === nothing && return nothing
    telemetry.exceptions += 1
    return nothing
end

@inline function record_scan_retry!(telemetry)
    telemetry = _require_solver_work_telemetry(telemetry)
    telemetry === nothing && return nothing
    telemetry.scan_retries += 1
    return nothing
end

function reset_solver_work!(telemetry::SolverWorkTelemetry)
    for field in fieldnames(SolverWorkTelemetry)
        setfield!(telemetry, field, 0)
    end
    return telemetry
end

function solver_work_snapshot(telemetry::SolverWorkTelemetry)
    return (; (field => getfield(telemetry, field) for field in fieldnames(SolverWorkTelemetry))...)
end
