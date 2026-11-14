"""Convergence evidence carried independently of the legacy acceptance flag.

`nothing` means the producing route did not record that item. `physical_valid`
checks finite state/thermodynamics, positive masses and the PNJL Polyakov range;
it does not assert global thermodynamic stability or replace custom hard rules.
"""
Base.@kwdef struct SolverConvergence
    solver_converged::Union{Nothing,Bool} = nothing
    residual_ok::Union{Nothing,Bool} = nothing
    physical_valid::Union{Nothing,Bool} = nothing
    residual_norm_max::Union{Nothing,Float64} = nothing
    source::Symbol = :not_recorded
    scan_promoted::Bool = false
end

@inline _result_convergence(result) = hasproperty(result, :convergence) ?
    _coerce_convergence(result.convergence) : SolverConvergence()
@inline _coerce_convergence(value::SolverConvergence) = value
@inline _coerce_convergence(::Nothing) = SolverConvergence()
function _coerce_convergence(value::NamedTuple)
    return SolverConvergence(; value...)
end

@inline function solver_convergence_view(result)
    c = result isa SolverConvergence ? result : _result_convergence(result)
    return (; c.solver_converged, c.residual_ok, c.physical_valid,
        c.residual_norm_max, c.source, c.scan_promoted)
end

function _convergence_physical_valid(result)
    return all(isfinite, result.solution) && all(isfinite, result.x_state) &&
        all(isfinite, result.mu_vec) && _thermo_quantities_finite(result) &&
        is_physical_solution(result.x_state, result.masses)
end

function _convergence_from_candidate(result; solver_converged=nothing,
    residual_norm_max::Real, source::Symbol)
    return SolverConvergence(;
        solver_converged,
        residual_ok=isfinite(result.residual_norm) && result.residual_norm <= residual_norm_max,
        physical_valid=_convergence_physical_valid(result),
        residual_norm_max=Float64(residual_norm_max), source,
    )
end

@inline function _record_gap_convergence!(destination, result)
    destination === nothing || (destination[] = Bool(result.f_converged))
    return nothing
end
