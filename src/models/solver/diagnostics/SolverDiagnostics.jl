const DIAGNOSTIC_LEVELS = (:none, :summary, :full)

@inline function _resolve_diagnostic_level(kwargs::Dict{Symbol,Any})::Symbol
    level = get(kwargs, :diagnostic_level, :none)
    level isa Symbol || throw(ArgumentError("diagnostic_level must be Symbol, got $(typeof(level))"))
    level in DIAGNOSTIC_LEVELS || throw(ArgumentError("invalid diagnostic_level: $(level), expected one of $(DIAGNOSTIC_LEVELS)"))
    return level
end

@inline function _diagnostic_attempt_origin(candidate)
    if hasproperty(candidate, :governed_attempt_origin)
        return Symbol(getproperty(candidate, :governed_attempt_origin))
    elseif hasproperty(candidate, :fixedrho_joint_attempt_origin)
        return Symbol(getproperty(candidate, :fixedrho_joint_attempt_origin))
    elseif hasproperty(candidate, :attempt_origin)
        return Symbol(getproperty(candidate, :attempt_origin))
    end
    return :fallback
end

@inline function _diagnostic_endpoint_cause(candidate)
    if hasproperty(candidate, :endpoint_cause)
        value = getproperty(candidate, :endpoint_cause)
        return value === nothing ? nothing : Symbol(value)
    end
    if !Bool(get(candidate, :converged, false))
        if hasproperty(candidate, :failed_constraints) && !isempty(getproperty(candidate, :failed_constraints))
            return :hard_constraint_failed
        end
        return :nonconvergence
    end
    return :converged
end

@inline function _solver_diagnostic_from_candidate(candidate; seed_source::Union{Symbol,Nothing})
    return to_namedtuple(_summary_from_candidate(candidate; seed_source=seed_source))
end

@inline function _attach_solver_diagnostic(
    result::NamedTuple,
    selected_candidate,
    candidates,
    diagnostic_level::Symbol;
    seed_source::Union{Symbol,Nothing}=nothing,
    selection_reason::Union{Nothing,Symbol}=nothing,
    selection_reason_source::Symbol=:problem_spec_selector,
)
    diagnostic_level === :none && return result
    summary_typed = _summary_from_candidate(
        selected_candidate;
        seed_source=seed_source,
        selection_reason=selection_reason,
        selection_reason_source=selection_reason_source,
    )
    summary = to_namedtuple(summary_typed)
    if diagnostic_level === :summary
        return merge(result, (diagnostic=summary,))
    end
    full_typed = SolverDiagnosticFull(
        summary_typed,
        [_summary_from_candidate(c; seed_source=seed_source, selection_reason_source=selection_reason_source) |> SolverDiagnosticCandidate for c in candidates],
    )
    return merge(result, (diagnostic=to_namedtuple(full_typed),))
end
