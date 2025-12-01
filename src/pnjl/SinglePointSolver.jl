module SinglePointSolver

using ..SeedCache: find_initial_seed, DEFAULT_SEED_PATH, seed_state
using ..Constants_PNJL: ħc_MeV_fm

export describe_solver, run_single_point

function describe_solver()
    return Dict(
        "id" => "pnjl-gap",
        "name" => "PNJL Gap (Anisotropic)",
        "description" => "单点各向异性 PNJL 能隙求解（依赖 coarse seed 暖启）",
        "params_schema" => Dict(
            "fields" => [
                Dict("name" => "T_mev", "type" => "number", "required" => true, "description" => "温度 MeV"),
                Dict("name" => "mu_mev", "type" => "number", "required" => false, "description" => "平均化学势 MeV"),
                Dict("name" => "rho", "type" => "number", "required" => false, "description" => "归一化密度 ρ/ρ₀"),
                Dict("name" => "xi", "type" => "number", "required" => false, "description" => "各向异性参数"),
            ],
        ),
    )
end

function _seed_summary(seed)
    return Dict(
        "T_MeV" => seed.T * ħc_MeV_fm,
        "mu_MeV" => seed.mu * ħc_MeV_fm,
        "rho" => seed.rho,
        "xi" => seed.xi,
        "phi" => Dict(
            "u" => seed.phi_u,
            "d" => seed.phi_d,
            "s" => seed.phi_s,
        ),
        "Phi" => Dict(
            "Phi1" => seed.Phi1,
            "Phi2" => seed.Phi2,
        ),
        "mu_components_MeV" => Dict(
            "mu_u" => seed.mu_u * ħc_MeV_fm,
            "mu_d" => seed.mu_d * ħc_MeV_fm,
            "mu_s" => seed.mu_s * ħc_MeV_fm,
        ),
    )
end

function run_single_point(params::Union{NamedTuple, AbstractDict, Dict}; seed_path::AbstractString = DEFAULT_SEED_PATH)
    nearest = find_initial_seed(params; path = seed_path)
    response = Dict(
        "status" => "not-implemented",
        "message" => "PNJL 求解器尚未接入，已返回最近邻初值供外部调用。",
        "distance" => nearest.distance,
        "used_axes" => nearest.used_axes,
        "seed" => _seed_summary(nearest.seed),
        "initial_state" => Dict(
            "x" => nearest.state.x,
            "mu" => nearest.state.mu,
            "phi" => nearest.state.phi,
            "polyakov" => nearest.state.Polyakov,
        ),
    )
    return response
end

end # module
