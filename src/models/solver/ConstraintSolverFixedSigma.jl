function _solve_constraint_fixedsigma(
    model::AbstractQCDModel,
    T_fm::Real,
    sigma_target::Real;
    seed_guess::AbstractVector,
    mu0::Real=0.2,
    solver_primary::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:forward, xtol=1e-9, ftol=1e-9),
    solver_secondary::Union{Nothing, AbstractGapSolver}=nothing,
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    residual_norm_max::Real=1e-6,
    nlsolve_method::Symbol=:trust_region,
    allow_legacy_fallback::Bool=false,
    rho0::Real,
    physicality_check::Function=((_, _) -> true),
    nlsolve_kwargs...,
)
    st_ref = Ref{Any}(nothing)
    x_state_ref = Ref{Any}(nothing)
    mu_vec_ref = Ref{Any}(nothing)
    pressure_ref = Ref{Any}(nothing)
    rho_norm_ref = Ref{Any}(nothing)
    entropy_ref = Ref{Any}(nothing)
    energy_ref = Ref{Any}(nothing)
    masses_ref = Ref{Any}(nothing)

    residual_fn! = (F, x) -> begin
        μ = x[1]
        μ_vec = normalize_mu_vec(μ)

        st = _solve_gap_with_outer_fallback(
            model,
            T_fm,
            μ;
            st_prev=st_ref[],
            seed_guess=seed_guess,
            solver_primary=solver_primary,
            solver_secondary=solver_secondary,
            residual_norm_max=residual_norm_max,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
        )

        if st === nothing
            _constraint_failure!(F)
            return nothing
        end

        x_state = _to_state_svec(st)
        pressure = -omega(model, x_state, T_fm, μ_vec; p_num=p_num, t_num=t_num, xi=xi)

        pressure_mu = μtrial -> -omega(model, x_state, T_fm, μtrial; p_num=p_num, t_num=t_num, xi=xi)
        rho_vec = ForwardDiff.gradient(pressure_mu, μ_vec)
        rho_norm = sum(rho_vec) / (3.0 * rho0)

        pressure_T = τ -> -omega(model, x_state, τ, μ_vec; p_num=p_num, t_num=t_num, xi=xi)
        entropy = ForwardDiff.derivative(pressure_T, T_fm)

        energy = -pressure + sum(μ_vec .* rho_vec) + T_fm * entropy
        masses = _mass_from_state(model, x_state)

        st_ref[] = st
        x_state_ref[] = x_state
        mu_vec_ref[] = _to_mu_svec(μ_vec)
        pressure_ref[] = pressure
        rho_norm_ref[] = rho_norm
        entropy_ref[] = entropy
        energy_ref[] = energy
        masses_ref[] = masses

        n_B = rho_norm * rho0
        if !isfinite(n_B) || abs(n_B) <= 1e-12
            _constraint_failure!(F)
            return nothing
        end

        F[1] = convert(eltype(F), entropy / n_B - sigma_target)
        return nothing
    end

    res = nlsolve(
        residual_fn!,
        [Float64(mu0)];
        autodiff=:forward,
        method=nlsolve_method,
        xtol=1e-9,
        ftol=1e-9,
        nlsolve_kwargs...,
    )

    _refresh_constraint_state!(residual_fn!, res.zero)

    gap_norm = _gap_norm_from_state(model, x_state_ref[], mu_vec_ref[], T_fm; xi=xi, p_num=p_num, t_num=t_num)
    n_B_ref = rho_norm_ref[] * rho0
    sigma_residual = if isfinite(n_B_ref) && abs(n_B_ref) > 1e-12
        abs(entropy_ref[] / n_B_ref - sigma_target)
    else
        Inf
    end
    residual_norm = max(gap_norm, sigma_residual)

    omega_val = -pressure_ref[]
    thermo_finite = isfinite(omega_val) && isfinite(pressure_ref[]) && isfinite(rho_norm_ref[]) && isfinite(entropy_ref[]) && isfinite(energy_ref[])
    phys = physicality_check(x_state_ref[], masses_ref[]) && thermo_finite
    converged = res.f_converged && phys && isfinite(residual_norm) && residual_norm <= max(Float64(residual_norm_max), 1e-3)

    _ = allow_legacy_fallback
    return (
        converged=converged,
        solution=_pack_solution(x_state_ref[], mu_vec_ref[]),
        x_state=x_state_ref[],
        mu_vec=mu_vec_ref[],
        omega=omega_val,
        pressure=pressure_ref[],
        rho_norm=rho_norm_ref[],
        entropy=entropy_ref[],
        energy=energy_ref[],
        masses=masses_ref[],
        iterations=res.iterations,
        residual_norm=residual_norm,
        legacy_fallback_used=false,
    )
end
