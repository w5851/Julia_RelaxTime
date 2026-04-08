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
    hard_ok = hasproperty(candidate, :hard_constraint_ok) ? Bool(getproperty(candidate, :hard_constraint_ok)) : nothing
    failed = hasproperty(candidate, :failed_constraints) ? Symbol.(getproperty(candidate, :failed_constraints)) : Symbol[]
    continuity = hasproperty(candidate, :continuity_distance) ? Float64(getproperty(candidate, :continuity_distance)) : nothing
    error_kind = hasproperty(candidate, :error_kind) ? Symbol(getproperty(candidate, :error_kind)) : :none
    error_msg = hasproperty(candidate, :error_msg) ? String(getproperty(candidate, :error_msg)) : ""
    selection_reason = hasproperty(candidate, :selection_reason) ? Symbol(getproperty(candidate, :selection_reason)) : :none
    return (
        attempt_origin=_diagnostic_attempt_origin(candidate),
        seed_source=seed_source,
        hard_constraint_ok=hard_ok,
        failed_constraints=failed,
        error_kind=error_kind,
        error_msg=error_msg,
        selection_reason=selection_reason,
        endpoint_cause=_diagnostic_endpoint_cause(candidate),
        continuity_distance=continuity,
    )
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
    summary_raw = _solver_diagnostic_from_candidate(selected_candidate; seed_source=seed_source)
    summary = if selection_reason === nothing
        (; summary_raw..., selection_reason_source=selection_reason_source)
    else
        (; summary_raw..., selection_reason=selection_reason, selection_reason_source=selection_reason_source)
    end
    if diagnostic_level === :summary
        return merge(result, (diagnostic=summary,))
    end
    full_candidates = [
        (; _solver_diagnostic_from_candidate(c; seed_source=seed_source)...,
            selection_reason_source=selection_reason_source,
        ) for c in candidates
    ]
    return merge(result, (diagnostic=(; summary..., candidates=full_candidates),))
end
