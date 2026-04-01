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
