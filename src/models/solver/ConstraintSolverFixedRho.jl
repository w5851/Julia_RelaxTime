function _solve_constraint_fixedrho(
    model::AbstractQCDModel,
    T_fm::Real,
    rho_target::Real;
    seed_guess::AbstractVector,
    mu0::Real=0.2,
    solver_primary::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:forward, xtol=1e-9, ftol=1e-9),
    solver_secondary::Union{Nothing, AbstractGapSolver}=nothing,
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    residual_norm_max::Real=1e-6,
    nlsolve_method::Symbol=:newton,
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

        F[1] = convert(eltype(F), rho_norm - rho_target)
        return nothing
    end

    function run_outer_attempt(mu_init::Real, method::Symbol)
        local res
        try
            res = _run_outer_nlsolve(
                residual_fn!,
                [mu_init];
                nlsolve_method=method,
                nlsolve_kwargs...,
            )
        catch err
            err isa InterruptException && rethrow()
            return nothing
        end

        if !_mode_outer_state_ready(st_ref, x_state_ref, mu_vec_ref, pressure_ref, rho_norm_ref, entropy_ref, energy_ref, masses_ref)
            return nothing
        end

        residual_norm = _compose_mode_residual_norm(
            model,
            x_state_ref[],
            mu_vec_ref[],
            T_fm,
            (rho_norm_ref[], rho_target);
            xi=xi,
            p_num=p_num,
            t_num=t_num,
        )

        omega_val = -pressure_ref[]
        thermo_finite = isfinite(omega_val) && isfinite(pressure_ref[]) && isfinite(rho_norm_ref[]) && isfinite(entropy_ref[]) && isfinite(energy_ref[])
        phys = physicality_check(x_state_ref[], masses_ref[]) && thermo_finite
        converged = res.f_converged && phys && isfinite(residual_norm) && residual_norm <= Float64(residual_norm_max)

        base = _build_mode_result_from_outer_state(
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
        )
        return (; base..., mode=method)
    end

    function evaluate_mu_candidate(mu_eval::Real)
        μ_vec = normalize_mu_vec(mu_eval)
        st = _solve_gap_with_outer_fallback(
            model,
            T_fm,
            mu_eval;
            st_prev=nothing,
            seed_guess=seed_guess,
            solver_primary=solver_primary,
            solver_secondary=solver_secondary,
            residual_norm_max=residual_norm_max,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
        )
        st === nothing && return nothing

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

        residual_norm = _compose_mode_residual_norm(
            model,
            x_state,
            μ_vec,
            T_fm,
            (thermo.rho_norm, rho_target);
            xi=xi,
            p_num=p_num,
            t_num=t_num,
        )

        thermo_finite = isfinite(thermo.omega) && isfinite(thermo.pressure) && isfinite(thermo.rho_norm) && isfinite(thermo.entropy) && isfinite(thermo.energy)
        phys = physicality_check(x_state, thermo.masses) && thermo_finite
        converged = phys && isfinite(residual_norm) && residual_norm <= Float64(residual_norm_max)

        base = _build_mode_result_from_outer_state(
            x_state,
            _to_mu_svec(μ_vec),
            thermo.pressure,
            thermo.rho_norm,
            thermo.entropy,
            thermo.energy,
            thermo.masses,
            residual_norm;
            iterations=0,
            converged=converged,
        )
        return (; base..., mode=:direct_mu)
    end

    @inline function prefer_candidate(cand, incumbent)
        if incumbent === nothing
            return true
        end
        if cand.converged != incumbent.converged
            return cand.converged
        end
        return cand.residual_norm < incumbent.residual_norm ||
               (cand.residual_norm == incumbent.residual_norm && cand.pressure > incumbent.pressure)
    end

    mu_seed = default_mu0_from_seed(seed_guess)
    attempt_specs = ((Float64(mu0), nlsolve_method),)

    best = nothing
    for (mu_init, method) in attempt_specs
        cand = run_outer_attempt(mu_init, method)
        cand === nothing && continue
        if prefer_candidate(cand, best)
            best = cand
        end
    end

    if best === nothing || !best.converged
        mu_min = 0.0
        mu_max = max(Float64(mu_seed), Float64(mu0), 2.0)
        mu_grid = collect(range(mu_min, mu_max; length=25))
        grid_candidates = NamedTuple[]
        for μ in mu_grid
            cand = evaluate_mu_candidate(μ)
            cand === nothing && continue
            push!(grid_candidates, cand)
        end

        if !isempty(grid_candidates)
            for cand in grid_candidates
                if prefer_candidate(cand, best)
                    best = cand
                end
            end

            sort!(grid_candidates; by=c -> c.mu_vec[1])
            for i in 1:(length(grid_candidates)-1)
                c1 = grid_candidates[i]
                c2 = grid_candidates[i + 1]
                f1 = c1.rho_norm - rho_target
                f2 = c2.rho_norm - rho_target
                if !(isfinite(f1) && isfinite(f2))
                    continue
                end
                if f1 == 0.0
                    best = c1
                    break
                end
                if f1 * f2 > 0.0
                    continue
                end

                left = c1.mu_vec[1]
                right = c2.mu_vec[1]
                f_left = f1
                for _ in 1:16
                    mid = 0.5 * (left + right)
                    cmid = evaluate_mu_candidate(mid)
                    cmid === nothing && break

                    if prefer_candidate(cmid, best)
                        best = cmid
                    end

                    f_mid = cmid.rho_norm - rho_target
                    if abs(f_mid) <= Float64(residual_norm_max)
                        break
                    end
                    if f_left * f_mid <= 0.0
                        right = mid
                    else
                        left = mid
                        f_left = f_mid
                    end
                end
                break
            end
        end
    end

    best === nothing && throw(ArgumentError("FixedRho outer solve failed for all attempts"))

    return (
        converged=best.converged,
        solution=best.solution,
        x_state=best.x_state,
        mu_vec=best.mu_vec,
        omega=best.omega,
        pressure=best.pressure,
        rho_norm=best.rho_norm,
        entropy=best.entropy,
        energy=best.energy,
        masses=best.masses,
        iterations=best.iterations,
        residual_norm=best.residual_norm,
    )
end
