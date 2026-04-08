function _solve_constraint_fixedmu(
    model::AbstractQCDModel,
    T_fm::Real,
    μ_fm::Real;
    seed_guess::AbstractVector,
    solver::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:forward, xtol=1e-9, ftol=1e-9),
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    residual_norm_max::Real=1e-6,
    physicality_check::Function=((_, _) -> true),
    seed_candidates::Union{Nothing, AbstractVector}=nothing,
    hard_constraints::Union{Nothing, AbstractVector{<:HardConstraintRule}}=nothing,
    allow_legacy_fallback::Bool=false,
    nlsolve_kwargs...,
)
    seed_pool = if seed_candidates === nothing
        _build_default_seed_candidates(seed_guess)
    else
        [Float64.(s) for s in seed_candidates]
    end

    rules = hard_constraints === nothing ? default_hard_constraint_rules(; physicality_check=physicality_check) : hard_constraints

    _ = allow_legacy_fallback
    _ = nlsolve_kwargs

    candidates = Main.Models.execute_attempt_pool(seed_pool;
        stop_on_first_success=true,
        evaluate_all_attempts=true,
        evaluate_attempt=(seed, seed_index) -> begin
            solver_pool = if solver isa NLsolveGapSolver && solver.method != :trust_region
                (
                    solver,
                    NLsolveGapSolver(
                        method=:trust_region,
                        xtol=solver.xtol,
                        ftol=solver.ftol,
                        jacobian=solver.jacobian,
                    ),
                )
            else
                (solver,)
            end

            for local_solver in solver_pool
                try
                    st = solve_gap(
                        model,
                        T_fm,
                        μ_fm;
                        solver_backend=:models,
                        solver=local_solver,
                        initial_guess=seed,
                        residual_norm_max=residual_norm_max,
                        xi=xi,
                        p_num=p_num,
                        t_num=t_num,
                    )

                    raw = _compute_fixedmu_candidate(
                        model,
                        T_fm,
                        μ_fm,
                        st,
                        residual_norm_max;
                        xi=xi,
                        p_num=p_num,
                        t_num=t_num,
                    )
                    ok, failed = evaluate_hard_constraints(raw, rules)
                    candidate = (; raw..., hard_constraint_ok=ok, failed_constraints=failed, converged=ok, seed_index=Int(seed_index))
                    return candidate, ok
                catch err
                    err isa InterruptException && rethrow()
                end
            end

            raw = _empty_candidate(; state_n=5, mu_n=3, solution_n=5, residual_norm_max=residual_norm_max)
            return (; raw..., seed_index=Int(seed_index)), false
        end,
        on_error=(_, seed_index, _) -> begin
            raw = _empty_candidate(; state_n=5, mu_n=3, solution_n=5, residual_norm_max=residual_norm_max)
            return (; raw..., seed_index=Int(seed_index)), false
        end,
    )
    selected = select_pressure_max_candidate(candidates)
    s = selected.selected_candidate

    return (
        converged=s.converged,
        solution=s.solution,
        x_state=s.x_state,
        mu_vec=s.mu_vec,
        omega=s.omega,
        pressure=s.pressure,
        rho_norm=s.rho_norm,
        entropy=s.entropy,
        energy=s.energy,
        masses=s.masses,
        iterations=s.iterations,
        residual_norm=s.residual_norm,
        hard_constraint_ok=s.hard_constraint_ok,
        failed_constraints=s.failed_constraints,
        selection_reason=selected.selection_reason,
        selected_index=selected.selected_index,
        candidate_count=length(candidates),
    )
end
