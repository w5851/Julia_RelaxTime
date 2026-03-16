function normalize_pm_seed_pair(seed_pair)
    hadron_seed0 = Float64.(collect(seed_pair.hadron_seed0))
    quark_seed0 = Float64.(collect(seed_pair.quark_seed0))
    continuity_mode = get(seed_pair, :continuity_mode, :branch_local)
    fallback_mode = get(seed_pair, :fallback_mode, :none)
    return PMSeedPair(
        hadron_seed0=hadron_seed0,
        quark_seed0=quark_seed0,
        continuity_mode=Symbol(continuity_mode),
        fallback_mode=Symbol(fallback_mode),
    )
end

function pm_next_seed_source(accepted::Bool, branch::Symbol)
    branch in (:hadron, :quark) || throw(ArgumentError("unsupported branch: $branch"))
    return accepted ? :previous_same_branch : :seed0
end

function _pm_solve_seed_equilibrium(T_MeV::Real, mu_MeV::Real, seed0::AbstractVector{<:Real};
        branch::Symbol=:hadron,
        xi::Real=0.0,
        solver_backend::Symbol=:legacy,
        p_num::Int=24,
        t_num::Int=8,
        residual_accept_tol::Float64=1e-6)
    T_fm = Float64(T_MeV) / Main.Constants_PNJL.ħc_MeV_fm
    mu_fm = Float64(mu_MeV) / Main.Constants_PNJL.ħc_MeV_fm

    if solver_backend === :models
        model = create_model(:PNJL)
        return solve_constraint(
            model,
            FixedMu(),
            T_fm;
            μ_fm=mu_fm,
            seed_guess=Float64.(seed0),
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            residual_norm_max=residual_accept_tol,
        )
    elseif solver_backend === :legacy
        phase_hint = branch === :quark ? :quark : :hadron
        strategy = DefaultSeed(Float64.(seed0), Float64.(seed0), phase_hint)
        return solve(
            FixedMu(),
            T_fm,
            mu_fm;
            xi=xi,
            seed_strategy=strategy,
            p_num=p_num,
            t_num=t_num,
        )
    end

    throw(ArgumentError("unsupported solver_backend: $solver_backend"))
end

function derive_pm_seed_pair(T_MeV::Real, mu_grid;
        xi::Real=0.0,
        solver_backend::Symbol=:legacy,
        p_num::Int=24,
        t_num::Int=8,
        residual_accept_tol::Float64=1e-6)
    mu_values = collect(Float64, mu_grid)
    isempty(mu_values) && throw(ArgumentError("mu_grid must not be empty"))

    hadron_seed = copy(HADRON_SEED_5)
    quark_seed = copy(QUARK_SEED_5)

    hadron_result = try
        _pm_solve_seed_equilibrium(T_MeV, first(mu_values), hadron_seed;
            branch=:hadron,
            xi=xi,
            solver_backend=solver_backend,
            p_num=p_num,
            t_num=t_num,
            residual_accept_tol=residual_accept_tol)
    catch
        nothing
    end

    if hadron_result !== nothing && Bool(getproperty(hadron_result, :converged))
        hadron_seed = Float64.(collect(getproperty(hadron_result, :x_state)))
    end

    quark_result = try
        _pm_solve_seed_equilibrium(T_MeV, last(mu_values), quark_seed;
            branch=:quark,
            xi=xi,
            solver_backend=solver_backend,
            p_num=p_num,
            t_num=t_num,
            residual_accept_tol=residual_accept_tol)
    catch
        nothing
    end

    if quark_result !== nothing && Bool(getproperty(quark_result, :converged))
        quark_seed = Float64.(collect(getproperty(quark_result, :x_state)))
    end

    return PMSeedPair(
        hadron_seed0=hadron_seed,
        quark_seed0=quark_seed,
        continuity_mode=:branch_local,
        fallback_mode=:none,
    )
end
