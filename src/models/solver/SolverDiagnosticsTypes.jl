"""SolverDiagnosticsTypes

诊断数据内部类型化表示，以及对旧 NamedTuple 视图的兼容映射。
"""

struct SolverDiagnosticSummary
    attempt_origin::Symbol
    seed_source::Union{Symbol,Nothing}
    hard_constraint_ok::Union{Bool,Nothing}
    failed_constraints::Vector{Symbol}
    error_kind::Symbol
    error_msg::String
    selection_reason::Symbol
    endpoint_cause::Union{Symbol,Nothing}
    continuity_distance::Union{Float64,Nothing}
    selection_reason_source::Symbol
end

struct SolverDiagnosticCandidate
    summary::SolverDiagnosticSummary
end

struct SolverDiagnosticFull
    summary::SolverDiagnosticSummary
    candidates::Vector{SolverDiagnosticCandidate}
end

@inline function to_namedtuple(summary::SolverDiagnosticSummary)
    return (
        attempt_origin=summary.attempt_origin,
        seed_source=summary.seed_source,
        hard_constraint_ok=summary.hard_constraint_ok,
        failed_constraints=summary.failed_constraints,
        error_kind=summary.error_kind,
        error_msg=summary.error_msg,
        selection_reason=summary.selection_reason,
        endpoint_cause=summary.endpoint_cause,
        continuity_distance=summary.continuity_distance,
        selection_reason_source=summary.selection_reason_source,
    )
end

@inline function to_namedtuple(candidate::SolverDiagnosticCandidate)
    return to_namedtuple(candidate.summary)
end

@inline function to_namedtuple(full::SolverDiagnosticFull)
    summary = to_namedtuple(full.summary)
    candidates = [to_namedtuple(c) for c in full.candidates]
    return (; summary..., candidates=candidates)
end

@inline function _summary_from_candidate(
    candidate;
    seed_source::Union{Symbol,Nothing}=nothing,
    selection_reason::Union{Nothing,Symbol}=nothing,
    selection_reason_source::Symbol=:problem_spec_selector,
)
    hard_ok = hasproperty(candidate, :hard_constraint_ok) ? Bool(getproperty(candidate, :hard_constraint_ok)) : nothing
    failed = hasproperty(candidate, :failed_constraints) ? Symbol.(getproperty(candidate, :failed_constraints)) : Symbol[]
    continuity = if hasproperty(candidate, :continuity_distance)
        value = getproperty(candidate, :continuity_distance)
        value === nothing ? nothing : Float64(value)
    else
        nothing
    end
    error_kind = hasproperty(candidate, :error_kind) ? Symbol(getproperty(candidate, :error_kind)) : :none
    error_msg = hasproperty(candidate, :error_msg) ? String(getproperty(candidate, :error_msg)) : ""
    reason = selection_reason === nothing ?
             (hasproperty(candidate, :selection_reason) ? Symbol(getproperty(candidate, :selection_reason)) : :none) :
             selection_reason
    endpoint = _diagnostic_endpoint_cause(candidate)
    return SolverDiagnosticSummary(
        _diagnostic_attempt_origin(candidate),
        seed_source,
        hard_ok,
        failed,
        error_kind,
        error_msg,
        reason,
        endpoint,
        continuity,
        selection_reason_source,
    )
end
