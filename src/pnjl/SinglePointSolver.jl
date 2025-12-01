module SinglePointSolver

using ..SeedCache: find_initial_seed, DEFAULT_SEED_PATH, seed_state
using ..Constants_PNJL: ħc_MeV_fm
using ..AnisoGapSolver: solve_fixed_rho, solve_fixed_mu, DEFAULT_RHO_GUESS, DEFAULT_MU_GUESS

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

_get_param(params, keys...) = begin
    for key in keys
        if params isa AbstractDict
            if haskey(params, key)
                return params[key]
            elseif key isa Symbol && haskey(params, Symbol(key))
                return params[Symbol(key)]
            end
        end
        if hasproperty(params, key)
            return getproperty(params, key)
        end
    end
    return nothing
end

function _normalize_request(params)
    T_mev = _get_param(params, :T_mev, :TMeV, :temperature_mev)
    T_mev === nothing && (T_mev = _get_param(params, :T, :temperature))
    T_mev === nothing && error("必须提供 T (MeV)")
    T_fm = T_mev / ħc_MeV_fm

    mu_val = _get_param(params, :mu_mev, :muMeV)
    if mu_val === nothing
        mu_nat = _get_param(params, :mu, :mu_fm)
        mu_val = mu_nat === nothing ? nothing : mu_nat * ħc_MeV_fm
    end

    rho = _get_param(params, :rho, :density)
    xi = _get_param(params, :xi, :anisotropy)
    xi = xi === nothing ? 0.0 : Float64(xi)

    mode = if rho !== nothing
        :rho
    elseif mu_val !== nothing
        :mu
    else
        error("需要提供 rho 或 mu 之一用来选择求解模式")
    end

    return Dict(
        :T_mev => Float64(T_mev),
        :T_fm => Float64(T_fm),
        :mu_mev => mu_val === nothing ? nothing : Float64(mu_val),
        :mu_fm => mu_val === nothing ? nothing : Float64(mu_val / ħc_MeV_fm),
        :rho => rho === nothing ? nothing : Float64(rho),
        :xi => xi,
        :mode => mode,
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

function _extract_seed_vectors(nearest, mode)
    if nearest === nothing
        return Float64.(DEFAULT_RHO_GUESS), Float64.(DEFAULT_MU_GUESS)
    end
    x_vec = nearest.state.x
    rho_seed = length(x_vec) == 8 ? Float64.(x_vec) : Float64.(DEFAULT_RHO_GUESS)
    mu_seed = length(x_vec) >= 5 ? Float64.(x_vec[1:5]) : Float64.(DEFAULT_MU_GUESS)
    if mode == :rho
        return rho_seed, mu_seed
    else
        return rho_seed, mu_seed
    end
end

function run_single_point(params::Union{NamedTuple, AbstractDict, Dict}; seed_path::AbstractString = DEFAULT_SEED_PATH, solver_kwargs...)
    request = _normalize_request(params)
    nearest = try
        find_initial_seed(params; path = seed_path)
    catch
        nothing
    end

    rho_seed, mu_seed = _extract_seed_vectors(nearest, request[:mode])
    xi = request[:xi]
    result = if request[:mode] == :rho
        solve_fixed_rho(request[:T_fm], request[:rho]; xi = xi, seed_state = rho_seed, solver_kwargs...)
    else
        solve_fixed_mu(request[:T_fm], request[:mu_fm]; xi = xi, seed_state = mu_seed, solver_kwargs...)
    end

    status = result.converged ? "ok" : "warning"
    seed_payload = nearest === nothing ? nothing : _seed_summary(nearest.seed)
    seed_distance = nearest === nothing ? nothing : nearest.distance
    mu_components = Dict(
        "u_fm" => result.mu_vec[1],
        "d_fm" => result.mu_vec[2],
        "s_fm" => result.mu_vec[3],
        "u_MeV" => result.mu_vec[1] * ħc_MeV_fm,
        "d_MeV" => result.mu_vec[2] * ħc_MeV_fm,
        "s_MeV" => result.mu_vec[3] * ħc_MeV_fm,
    )

    payload = Dict(
        "status" => status,
        "mode" => String(result.mode),
        "xi" => result.xi,
        "converged" => result.converged,
        "iterations" => result.iterations,
        "residual_norm" => result.residual_norm,
        "thermo" => Dict(
            "pressure_fm4" => result.pressure,
            "rho" => result.rho,
            "entropy_fm3" => result.entropy,
            "energy_fm4" => result.energy,
        ),
        "solution_vector" => result.solution,
        "mu_components" => mu_components,
        "seed" => seed_payload,
        "seed_distance" => seed_distance,
    )
    return payload
end

end # module
