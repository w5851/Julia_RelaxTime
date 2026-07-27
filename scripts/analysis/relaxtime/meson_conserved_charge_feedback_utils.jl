"""
Pure numerical helpers for the charged-meson conserved-charge feedback spike.

The module is analysis-only. It does not solve a PNJL state or evaluate a BU
kernel itself; callers inject one candidate evaluator so the outer algorithm
and charge/sign contracts remain unit-testable.
"""
module MesonConservedChargeFeedbackUtils

using LinearAlgebra: norm, svdvals

export charged_meson_chemical_potentials
export charged_meson_conserved_densities
export total_conserved_charge_residual
export solve_outer_conserved_charge_feedback
export choose_freezeout_sqrts_grid
export feedback_initial_mu

@inline function _finite(name::AbstractString, value::Real)::Float64
    x = Float64(value)
    isfinite(x) || throw(ArgumentError("$(name) must be finite, got $(value)"))
    return x
end

@inline function _nonnegative(name::AbstractString, value::Real)::Float64
    x = _finite(name, value)
    x >= 0.0 || throw(ArgumentError("$(name) must be nonnegative, got $(value)"))
    return x
end

@inline function _required_real(payload, name::Symbol)::Float64
    hasproperty(payload, name) || throw(ArgumentError("candidate payload missing $(name)"))
    return _finite(String(name), getproperty(payload, name))
end

"""
Choose the largest diagnostic freeze-out grid that fits an approximate wall-time
budget.  The returned values are always ascending in `sqrt(s_NN)` and are kept
deliberately sparse because this helper is not a production scan planner.
"""
function choose_freezeout_sqrts_grid(
    median_seconds::Real;
    budget_seconds::Real=600.0,
)
    median_s = _finite("median_seconds", median_seconds)
    median_s > 0.0 || throw(ArgumentError("median_seconds must be positive"))
    budget_s = _finite("budget_seconds", budget_seconds)
    budget_s > 0.0 || throw(ArgumentError("budget_seconds must be positive"))

    grids = (
        Float64[3.0, 4.5, 7.7, 11.5, 20.0, 62.4, 200.0],
        Float64[3.0, 7.7, 11.5, 39.0, 200.0],
        Float64[3.0, 7.7, 200.0],
    )
    for grid in grids
        median_s * length(grid) <= budget_s && return copy(grid)
    end
    return copy(last(grids))
end

"""Select continuation `(mu_Q,mu_S)` or fall back to the current quark-only seed."""
function feedback_initial_mu(
    previous_mu_Q::Real,
    previous_mu_S::Real,
    previous_converged::Bool,
    fallback_mu_Q::Real,
    fallback_mu_S::Real,
)
    previous_converged && return (
        mu_Q=_finite("previous_mu_Q", previous_mu_Q),
        mu_S=_finite("previous_mu_S", previous_mu_S),
        source=:previous_feedback,
    )
    return (
        mu_Q=_finite("fallback_mu_Q", fallback_mu_Q),
        mu_S=_finite("fallback_mu_S", fallback_mu_S),
        source=:quark_only_fallback,
    )
end

"""Map `(mu_Q,mu_S)` to chemical-equilibrium charged-meson potentials."""
function charged_meson_chemical_potentials(mu_Q::Real, mu_S::Real)
    q = _finite("mu_Q", mu_Q)
    s = _finite("mu_S", mu_S)
    return (
        mu_pi_plus=q,
        mu_pi_minus=-q,
        mu_K_plus=q + s,
        mu_K_minus=-(q + s),
    )
end

"""Map four charged-meson number densities to net electric/strange density."""
function charged_meson_conserved_densities(
    n_pi_plus::Real,
    n_pi_minus::Real,
    n_K_plus::Real,
    n_K_minus::Real,
)
    pi_plus = _nonnegative("n_pi_plus", n_pi_plus)
    pi_minus = _nonnegative("n_pi_minus", n_pi_minus)
    K_plus = _nonnegative("n_K_plus", n_K_plus)
    K_minus = _nonnegative("n_K_minus", n_K_minus)
    return (
        rho_B=0.0,
        rho_Q=pi_plus - pi_minus + K_plus - K_minus,
        rho_S=K_plus - K_minus,
    )
end

function charged_meson_conserved_densities(payload)
    return charged_meson_conserved_densities(
        _required_real(payload, :n_pi_plus),
        _required_real(payload, :n_pi_minus),
        _required_real(payload, :n_K_plus),
        _required_real(payload, :n_K_minus),
    )
end

"""
Build affine quark+meson charge residuals, normalized by `rho0`.

The charge equation is never evaluated as `rho_Q/rho_B`; this keeps the
residual defined when the baryon density is small.
"""
function total_conserved_charge_residual(
    quark_densities,
    meson_densities;
    charge_to_baryon_ratio::Real=0.4,
    strangeness_density_target::Real=0.0,
    rho0::Real=0.16,
)
    rho_B_q = _required_real(quark_densities, :rho_B)
    rho_Q_q = _required_real(quark_densities, :rho_Q)
    rho_S_q = _required_real(quark_densities, :rho_S)
    rho_Q_M = _required_real(meson_densities, :rho_Q)
    rho_S_M = _required_real(meson_densities, :rho_S)
    target = _finite("charge_to_baryon_ratio", charge_to_baryon_ratio)
    rho_S_target = _finite("strangeness_density_target", strangeness_density_target)
    rho0_f = _finite("rho0", rho0)
    rho0_f > 0.0 || throw(ArgumentError("rho0 must be positive, got $(rho0)"))

    charge_raw = rho_Q_q + rho_Q_M - target * rho_B_q
    strangeness_raw = rho_S_q + rho_S_M - rho_S_target
    residual = Float64[charge_raw / rho0_f, strangeness_raw / rho0_f]
    return (
        charge_raw=charge_raw,
        strangeness_raw=strangeness_raw,
        normalized=residual,
        norm=norm(residual),
        rho_B_total=rho_B_q,
        rho_Q_total=rho_Q_q + rho_Q_M,
        rho_S_total=rho_S_q + rho_S_M,
    )
end

@inline function _candidate_residual(
    evaluator::Function,
    x::Vector{Float64};
    charge_to_baryon_ratio::Float64,
    strangeness_density_target::Float64,
    rho0::Float64,
)
    payload = evaluator(x[1], x[2])
    quark = (
        rho_B=_required_real(payload, :rho_B_q),
        rho_Q=_required_real(payload, :rho_Q_q),
        rho_S=_required_real(payload, :rho_S_q),
    )
    meson = charged_meson_conserved_densities(payload)
    residual = total_conserved_charge_residual(
        quark,
        meson;
        charge_to_baryon_ratio=charge_to_baryon_ratio,
        strangeness_density_target=strangeness_density_target,
        rho0=rho0,
    )
    return payload, meson, residual
end

@inline function _failure_result(
    reason::Symbol,
    x::Vector{Float64},
    evaluations::Int,
    iterations::Int,
    history,
    payload,
    meson,
    residual;
    message::AbstractString="",
)
    return (
        converged=false,
        reason=reason,
        message=String(message),
        mu_Q=x[1],
        mu_S=x[2],
        residual=residual,
        residual_norm=residual === nothing ? Inf : residual.norm,
        iterations=iterations,
        evaluation_count=evaluations,
        payload=payload,
        meson_densities=meson,
        history=history,
    )
end

"""
    solve_outer_conserved_charge_feedback(evaluator, initial_mu_Q, initial_mu_S; ...)

Solve the two partial-feedback charge equations with a damped finite-difference
Newton method. `evaluator(mu_Q,mu_S)` must return quark B/Q/S densities plus
`n_pi_plus`, `n_pi_minus`, `n_K_plus`, and `n_K_minus`.
"""
function solve_outer_conserved_charge_feedback(
    evaluator::Function,
    initial_mu_Q::Real,
    initial_mu_S::Real;
    charge_to_baryon_ratio::Real=0.4,
    strangeness_density_target::Real=0.0,
    rho0::Real=0.16,
    residual_tolerance::Real=1e-4,
    finite_difference_step::Real=1e-3,
    maximum_step::Real=0.25,
    step_tolerance::Real=1e-8,
    jacobian_singular_tolerance::Real=1e-10,
    max_iterations::Int=8,
    max_evaluations::Int=40,
    line_search_factors::Tuple=(1.0, 0.5, 0.25, 0.125),
)
    x = Float64[_finite("initial_mu_Q", initial_mu_Q), _finite("initial_mu_S", initial_mu_S)]
    target = _finite("charge_to_baryon_ratio", charge_to_baryon_ratio)
    rho_S_target = _finite("strangeness_density_target", strangeness_density_target)
    rho0_f = _finite("rho0", rho0)
    rho0_f > 0.0 || throw(ArgumentError("rho0 must be positive, got $(rho0)"))
    tolerance = _finite("residual_tolerance", residual_tolerance)
    tolerance > 0.0 || throw(ArgumentError("residual_tolerance must be positive"))
    fd_step = _finite("finite_difference_step", finite_difference_step)
    fd_step > 0.0 || throw(ArgumentError("finite_difference_step must be positive"))
    max_step = _finite("maximum_step", maximum_step)
    max_step > 0.0 || throw(ArgumentError("maximum_step must be positive"))
    step_tol = _finite("step_tolerance", step_tolerance)
    step_tol >= 0.0 || throw(ArgumentError("step_tolerance must be nonnegative"))
    singular_tol = _finite("jacobian_singular_tolerance", jacobian_singular_tolerance)
    singular_tol > 0.0 || throw(ArgumentError("jacobian_singular_tolerance must be positive"))
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    max_evaluations >= 3 || throw(ArgumentError("max_evaluations must be at least 3"))
    isempty(line_search_factors) && throw(ArgumentError("line_search_factors cannot be empty"))
    factors = Float64.([line_search_factors...])
    all(a -> isfinite(a) && 0.0 < a <= 1.0, factors) || throw(ArgumentError(
        "line_search_factors must lie in (0,1]",
    ))

    history = NamedTuple[]
    evaluations = 0
    payload = nothing
    meson = nothing
    residual = nothing

    try
        payload, meson, residual = _candidate_residual(
            evaluator,
            x;
            charge_to_baryon_ratio=target,
            strangeness_density_target=rho_S_target,
            rho0=rho0_f,
        )
        evaluations += 1
    catch err
        err isa InterruptException && rethrow()
        return _failure_result(
            :initial_evaluation_failed,
            x,
            evaluations + 1,
            0,
            history,
            nothing,
            nothing,
            nothing;
            message=sprint(showerror, err),
        )
    end
    push!(history, (iteration=0, mu_Q=x[1], mu_S=x[2], residual_norm=residual.norm))

    for iteration in 1:max_iterations
        if residual.norm <= tolerance
            return (
                converged=true,
                reason=:residual_tolerance,
                message="",
                mu_Q=x[1],
                mu_S=x[2],
                residual=residual,
                residual_norm=residual.norm,
                iterations=iteration - 1,
                evaluation_count=evaluations,
                payload=payload,
                meson_densities=meson,
                history=history,
            )
        end
        evaluations + 2 <= max_evaluations || return _failure_result(
            :evaluation_budget,
            x,
            evaluations,
            iteration - 1,
            history,
            payload,
            meson,
            residual,
        )

        J = zeros(2, 2)
        for j in 1:2
            trial = copy(x)
            h = fd_step * max(1.0, abs(x[j]))
            trial[j] += h
            try
                _, _, shifted = _candidate_residual(
                    evaluator,
                    trial;
                    charge_to_baryon_ratio=target,
                    strangeness_density_target=rho_S_target,
                    rho0=rho0_f,
                )
                evaluations += 1
                J[:, j] .= (shifted.normalized .- residual.normalized) ./ h
            catch err
                err isa InterruptException && rethrow()
                return _failure_result(
                    :jacobian_evaluation_failed,
                    x,
                    evaluations + 1,
                    iteration - 1,
                    history,
                    payload,
                    meson,
                    residual;
                    message=sprint(showerror, err),
                )
            end
        end

        singular_values = svdvals(J)
        if !all(isfinite, singular_values) || minimum(singular_values) <= singular_tol
            return _failure_result(
                :singular_jacobian,
                x,
                evaluations,
                iteration - 1,
                history,
                payload,
                meson,
                residual,
            )
        end

        step = -(J \ residual.normalized)
        step_norm = norm(step)
        if step_norm > max_step
            step .*= max_step / step_norm
            step_norm = max_step
        end

        accepted = false
        for alpha in factors
            evaluations < max_evaluations || break
            trial = x .+ alpha .* step
            try
                trial_payload, trial_meson, trial_residual = _candidate_residual(
                    evaluator,
                    trial;
                    charge_to_baryon_ratio=target,
                    strangeness_density_target=rho_S_target,
                    rho0=rho0_f,
                )
                evaluations += 1
                if trial_residual.norm < residual.norm
                    x = trial
                    payload = trial_payload
                    meson = trial_meson
                    residual = trial_residual
                    accepted = true
                    push!(history, (
                        iteration=iteration,
                        mu_Q=x[1],
                        mu_S=x[2],
                        residual_norm=residual.norm,
                    ))
                    alpha * step_norm <= step_tol && residual.norm > tolerance && return _failure_result(
                        :step_stalled,
                        x,
                        evaluations,
                        iteration,
                        history,
                        payload,
                        meson,
                        residual,
                    )
                    break
                end
            catch err
                err isa InterruptException && rethrow()
                evaluations += 1
            end
        end

        accepted || return _failure_result(
            evaluations >= max_evaluations ? :evaluation_budget : :line_search_failed,
            x,
            evaluations,
            iteration - 1,
            history,
            payload,
            meson,
            residual,
        )
    end

    return _failure_result(
        :max_iterations,
        x,
        evaluations,
        max_iterations,
        history,
        payload,
        meson,
        residual,
    )
end

end # module MesonConservedChargeFeedbackUtils
