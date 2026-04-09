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
        thermo = _compute_mode_thermo_quantities(
            model,
            x_state,
            T_fm,
            μ_vec;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            rho0_scale=rho0,
        )

        st_ref[] = st
        x_state_ref[] = x_state
        mu_vec_ref[] = _to_mu_svec(μ_vec)
        pressure_ref[] = thermo.pressure
        rho_norm_ref[] = thermo.rho_norm
        entropy_ref[] = thermo.entropy
        energy_ref[] = thermo.energy
        masses_ref[] = thermo.masses

        n_B = thermo.rho_norm * rho0
        if !isfinite(n_B) || abs(n_B) <= 1e-12
            _constraint_failure!(F)
            return nothing
        end

        F[1] = convert(eltype(F), thermo.entropy / n_B - sigma_target)
        return nothing
    end

    res = _run_outer_nlsolve(
        residual_fn!,
        [mu0];
        nlsolve_method=nlsolve_method,
        nlsolve_kwargs...,
    )

    n_B_ref = rho_norm_ref[] * rho0
    sigma_ratio = if isfinite(n_B_ref) && abs(n_B_ref) > 1e-12
        entropy_ref[] / n_B_ref
    else
        Inf
    end
    residual_norm = _compose_mode_residual_norm(
        model,
        x_state_ref[],
        mu_vec_ref[],
        T_fm,
        (sigma_ratio, sigma_target);
        xi=xi,
        p_num=p_num,
        t_num=t_num,
    )

    omega_val = -pressure_ref[]
    thermo_finite = isfinite(omega_val) && isfinite(pressure_ref[]) && isfinite(rho_norm_ref[]) && isfinite(entropy_ref[]) && isfinite(energy_ref[])
    phys = physicality_check(x_state_ref[], masses_ref[]) && thermo_finite
    converged = res.f_converged && phys && isfinite(residual_norm) && residual_norm <= max(Float64(residual_norm_max), 1e-3)

    return _build_mode_result_from_outer_state(
        x_state_ref[],
        mu_vec_ref[],
        pressure_ref[],
        rho_norm_ref[],
        entropy_ref[],
        energy_ref[],
        masses_ref[],
        residual_norm;
        iterations=res.iterations,
        converged=converged,
        legacy_fallback_used=false,
    )
end
