# Shared equilibrium stage of the trho_asymmetric scan. The caller supplies
# Models; the PALC diagnostic reuses this stage without running meson densities.
function _trho_asymmetric_seed_candidates(T_fm::Float64, mode, equilibrium_seed)
    multiseeds = Models.get_all_seeds(Models.MultiSeed(), [T_fm], mode)
    seed_pool = if equilibrium_seed === nothing
        Models.build_seed_pool(mode;
            primary_seed=multiseeds[1],
            extra_seed_pool=multiseeds[2:end],
            seed_extend=(seed, _) -> Float64.(seed),
        )
    else
        Models.build_seed_pool(mode;
            primary_seed=equilibrium_seed,
            extra_seed_pool=multiseeds,
            seed_extend=(seed, _) -> Float64.(seed),
        )
    end
    return [entry.seed for entry in seed_pool]
end

function _solve_trho_asymmetric_equilibrium(model, mode, T_fm;
    seeds, xi, p_num, t_num, iterations, kwargs...)
    return Models.solve_multi(model, mode, T_fm;
        seeds=seeds, xi=xi, p_num=p_num, t_num=t_num, iterations=iterations,
        semantic_mode=:ground_state, evaluate_all_attempts=true, kwargs...)
end
