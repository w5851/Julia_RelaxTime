"""Three-field PNJL gap solve with the exact zero-temperature grand potential."""

if !isdefined(@__MODULE__, :ValidationUtils)
    include(joinpath(@__DIR__, "..", "..", "..", "utils", "ValidationUtils.jl"))
end
using .ValidationUtils: require_finite, require_positive_finite,
    require_nonnegative_finite, require_positive_integer
using LinearAlgebra: Symmetric, eigmin, opnorm

struct ZeroTemperatureGapResidual{F}
    potential::F
end

@inline function (residual::ZeroTemperatureGapResidual)(out, phi)
    gradient = ForwardDiff.gradient(residual.potential, SVector{3}(phi))
    copyto!(out, gradient)
    return nothing
end

function _zero_temperature_gap_potential(model::PNJLModel, mu_vec;
    xi=0.0, p_num::Int=64, t_num::Int=8, kwargs...)
    controls = (; thermo_quadrature_policy=:rs_reduced_adaptive, kwargs...)
    return GapOmegaFn{3, typeof(model), Float64, typeof(mu_vec), typeof(controls), typeof(xi)}(
        model, 0.0, mu_vec, p_num, t_num, xi, controls)
end

function _zero_temperature_seed(seed)
    if seed isa AbstractVector
        length(seed) in (3, 5) || throw(ArgumentError(
            "zero-temperature initial_guess must have 3 chiral or 5 mean-field components"))
    end
    phi = _as_seed_vec3(seed)
    foreach(value -> require_finite("initial_guess chiral component", value), phi)
    return phi
end

function _zero_temperature_seeds()
    # Cold dense matter keeps a large strange condensate; the high-temperature
    # guesses alone can converge to the intervening saddle instead of that root.
    seeds = [_zero_temperature_seed(SeedStrategies.HADRON_SEED_5),
        _zero_temperature_seed(SeedStrategies.HIGH_DENSITY_SEED_5)]
    append!(seeds, _zero_temperature_seed.(_pnjl_multi_seeds()))
    return unique(seeds)
end

function _zero_temperature_nlsolve(residual!, seed, solver, residual_norm_max, iterations)
    result = nlsolve(residual!, copy(seed); autodiff=solver.jacobian,
        method=solver.method, xtol=solver.xtol, ftol=solver.ftol, iterations)
    resid = Float64(result.residual_norm)
    result.f_converged && isfinite(resid) && resid <= residual_norm_max ||
        error("zero-temperature solve_gap did not converge (residual_norm=$resid, f_converged=$(result.f_converged))")
    return MeanFieldState(SVector{3, Float64}(result.zero); Phi=0.0, PhiBar=0.0)
end

function _try_zero_temperature_nlsolve(residual!, seed, solver, residual_norm_max, iterations)
    try
        return _zero_temperature_nlsolve(residual!, seed, solver, residual_norm_max, iterations)
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
end

function _zero_temperature_local_minimum(model, state, potential)
    all(isfinite, state.phi) || return false
    masses = calculate_mass_vec(model, state.phi)
    all(m -> isfinite(m) && m > 0, masses) || return false
    hessian = ForwardDiff.hessian(potential, state.phi)
    all(isfinite, hessian) || return false
    # Permit only roundoff-scale negative curvature when selecting a minimum.
    roundoff = 100eps(Float64) * max(1.0, opnorm(hessian, Inf))
    return eigmin(Symmetric(hessian)) >= -roundoff
end

"""Solve the chiral subsystem; explicit seeds retain stationary-branch semantics.

Without a seed, reuse the existing catalog's distinct chiral guesses and select
the smallest grand potential among converged physical local minima. Polyakov
variables are absent from this nonlinear system and returned as zeros.
"""
function _solve_gap_zero_temperature(model::PNJLModel, mu_vec;
    solver::AbstractGapSolver=NLsolveGapSolver(),
    initial_guess=nothing,
    residual_norm_max::Real=1e-6,
    iterations::Int=1000,
    xi::Real=0.0,
    p_num::Int=PNJLCore.DEFAULT_MOMENTUM_COUNT,
    t_num::Int=PNJLCore.DEFAULT_THETA_COUNT,
    thermo_quadrature_policy::Symbol=:rs_reduced_adaptive,
    thermo_quadrature_rtol::Float64=1e-8,
    thermo_quadrature_atol::Float64=1e-10,
    thermo_quadrature_maxevals::Int=10^7,
    kwargs...)
    mu = SVector{3, Float64}(normalize_mu_vec(mu_vec))
    foreach(value -> require_finite("mu_vec component", value), mu)
    require_positive_integer("iterations", iterations)
    require_positive_integer("p_num", p_num)
    require_positive_integer("t_num", t_num)
    require_positive_finite("residual_norm_max", residual_norm_max)
    PNJLIntegrals.validate_rs_anisotropy(xi)
    PNJLIntegrals.validate_thermal_quadrature_policy(thermo_quadrature_policy)
    thermo_quadrature_policy === :rs_reduced_adaptive || throw(ArgumentError(
        "zero-temperature solve_gap uses the analytic Fermi sea; select :rs_reduced_adaptive"))
    PNJLIntegrals.validate_thermal_quadrature_controls(
        thermo_quadrature_rtol, thermo_quadrature_atol, thermo_quadrature_maxevals)
    s = solver isa NLsolveGapSolver ? solver : NLsolveGapSolver()
    s.jacobian in (:forward, :finite) || throw(ArgumentError(
        "solver.jacobian must be :forward or :finite"))
    require_nonnegative_finite("solver.xtol", s.xtol)
    require_positive_finite("solver.ftol", s.ftol)
    potential = _zero_temperature_gap_potential(model, mu;
        xi, p_num, t_num, thermo_quadrature_rtol, thermo_quadrature_atol,
        thermo_quadrature_maxevals, kwargs...)
    residual! = ZeroTemperatureGapResidual(potential)

    if initial_guess !== nothing
        return _zero_temperature_nlsolve(residual!, _zero_temperature_seed(initial_guess),
            s, Float64(residual_norm_max), iterations)
    end

    seeds = _zero_temperature_seeds()
    fallback = NLsolveGapSolver(method=:trust_region, jacobian=s.jacobian, xtol=s.xtol, ftol=s.ftol)
    best_state = nothing
    best_omega = Inf
    for seed in seeds
        state = _try_zero_temperature_nlsolve(residual!, seed, s, residual_norm_max, iterations)
        if (state === nothing || !_zero_temperature_local_minimum(model, state, potential)) &&
            s.method !== :trust_region
            state = _try_zero_temperature_nlsolve(residual!, seed, fallback, residual_norm_max, iterations)
        end
        state === nothing && continue
        _zero_temperature_local_minimum(model, state, potential) || continue
        value = potential(state.phi)
        if isfinite(value) && value < best_omega
            best_state, best_omega = state, value
        end
    end
    best_state === nothing && error(
        "zero-temperature solve_gap found no physical local minimum from $(length(seeds)) chiral seeds")
    return best_state
end
