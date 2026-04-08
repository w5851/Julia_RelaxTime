function _solve_constraint_fixedasymrho(
    model::AbstractQCDModel,
    T_fm::Real,
    rho_target::Real,
    ud_ratio_target::Real,
    s_target::Real;
    seed_guess::AbstractVector,
    mu0::AbstractVector=Float64[0.2, 0.2, 0.2],
    solver_primary::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:forward, xtol=1e-9, ftol=1e-9),
    solver_secondary::Union{Nothing, AbstractGapSolver}=nothing,
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    residual_norm_max::Real=1e-6,
    nlsolve_method::Symbol=:trust_region,
    allow_legacy_fallback::Bool=false,
    rho0::Real,
    enforce_physicality::Bool=false,
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
    rho_vec_ref = Ref{Any}(nothing)

    residual_fn! = (F, x) -> begin
        μ_vec = SVector{3}(x[1], x[2], x[3])

        st = _solve_gap_with_outer_fallback(
            model,
            T_fm,
            μ_vec;
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

        rho_u, rho_d, rho_s = rho_vec[1], rho_vec[2], rho_vec[3]
        nB = sum(rho_vec) / (3.0 * rho0)
        ud_ratio = if abs(rho_d) > 1e-12
            rho_u / rho_d
        else
            rho_u / (rho_d >= 0 ? 1e-12 : -1e-12)
        end

        st_ref[] = st
        x_state_ref[] = x_state
        mu_vec_ref[] = μ_vec
        pressure_ref[] = pressure
        rho_norm_ref[] = rho_norm
        entropy_ref[] = entropy
        energy_ref[] = energy
        masses_ref[] = masses
        rho_vec_ref[] = _to_mu_svec(rho_vec)

        F[1] = convert(eltype(F), nB - rho_target)
        F[2] = convert(eltype(F), ud_ratio - ud_ratio_target)
        F[3] = convert(eltype(F), rho_s - s_target)
        return nothing
    end

    res = nlsolve(
        residual_fn!,
        Float64.(mu0);
        autodiff=:forward,
        method=nlsolve_method,
        xtol=1e-9,
        ftol=1e-9,
        nlsolve_kwargs...,
    )

    _refresh_constraint_state!(residual_fn!, res.zero)

    gap_norm = _gap_norm_from_state(model, x_state_ref[], mu_vec_ref[], T_fm; xi=xi, p_num=p_num, t_num=t_num)
    rho_u, rho_d, rho_s = rho_vec_ref[][1], rho_vec_ref[][2], rho_vec_ref[][3]
    nB = sum(rho_vec_ref[]) / (3.0 * rho0)
    ud_ratio = if abs(rho_d) > 1e-12
        rho_u / rho_d
    else
        rho_u / (rho_d >= 0 ? 1e-12 : -1e-12)
    end
    nB_residual = abs(nB - rho_target)
    ud_residual = abs(ud_ratio - ud_ratio_target)
    s_residual = abs(rho_s - s_target)
    residual_norm = max(gap_norm, nB_residual, ud_residual, s_residual)

    omega_val = -pressure_ref[]
    thermo_finite = isfinite(omega_val) && isfinite(pressure_ref[]) && isfinite(rho_norm_ref[]) && isfinite(entropy_ref[]) && isfinite(energy_ref[])
    phys = physicality_check(x_state_ref[], masses_ref[]) && thermo_finite
    converged = if enforce_physicality
        res.f_converged && phys && isfinite(residual_norm) && residual_norm <= max(Float64(residual_norm_max), 1e-3)
    else
        isfinite(residual_norm) && residual_norm <= max(Float64(residual_norm_max), 1e-3)
    end

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
