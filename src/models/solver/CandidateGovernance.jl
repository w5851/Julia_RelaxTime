"""
    Candidate governance context and adapters.
"""

const HardRule = Function
const CandidateSelector = Function

@inline function build_candidate_context(
    mode::ConstraintMode;
    continuity_seed_available::Bool=false,
    phase_hint::Symbol=:unknown,
    residual_norm_max::Real=1e-6,
    prefer_continuity::Bool=false,
)
    return (
        mode=mode,
        continuity_seed_available=Bool(continuity_seed_available),
        phase_hint=phase_hint,
        residual_norm_max=Float64(residual_norm_max),
        prefer_continuity=Bool(prefer_continuity),
    )
end

function execute_attempt_pool(
    attempts;
    evaluate_attempt::Function,
    on_error::Function,
    stop_on_first_success::Bool=true,
    evaluate_all_attempts::Bool=false,
)
    candidates = NamedTuple[]
    for (attempt_index, attempt_cfg) in enumerate(attempts)
        candidate, success = try
            evaluate_attempt(attempt_cfg, attempt_index)
        catch err
            err isa InterruptException && rethrow()
            on_error(attempt_cfg, attempt_index, err)
        end
        push!(candidates, candidate)
        if !evaluate_all_attempts && stop_on_first_success && Bool(success)
            break
        end
    end
    return candidates
end
