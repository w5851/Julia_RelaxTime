module SeedCache

using DelimitedFiles
using ..Constants_PNJL: ħc_MeV_fm

export SeedPoint, load_seed_table, find_initial_seed, seed_state, DEFAULT_SEED_PATH

struct SeedPoint
    T::Float64        # 温度 (fm⁻¹)
    mu::Float64       # 公共化学势 (fm⁻¹)
    rho::Float64      # 归一化密度 (单位 rho0)
    xi::Float64       # 各向异性参数
    phi_u::Float64
    phi_d::Float64
    phi_s::Float64
    Phi1::Float64
    Phi2::Float64
    mu_u::Float64
    mu_d::Float64
    mu_s::Float64
    source::String
end

const DEFAULT_SEED_PATH = normpath(joinpath(@__DIR__, "..", "..", "data", "raw", "pnjl", "seeds", "aniso_coarse.csv"))
const DEFAULT_WEIGHTS = (T = 1.0, mu = 1.0, rho = 0.6, xi = 0.4)

const _seed_cache = Ref{Vector{SeedPoint}}(nothing)
const _seed_cache_path = Ref{String}("")

# ---------------- helper functions ----------------

_get_param(params, key, default) =
    params isa AbstractDict ? get(params, key, default) :
    hasproperty(params, key) ? getproperty(params, key) :
    default

function _normalize_request(params)
    T = let val = _get_param(params, :T_fm, nothing)
        val === nothing ? _get_param(params, :T, nothing) : val
    end
    if T === nothing
        T_mev = _get_param(params, :T_mev, nothing)
        T = T_mev === nothing ? nothing : T_mev / ħc_MeV_fm
    end

    mu = let val = _get_param(params, :mu_fm, nothing)
        val === nothing ? _get_param(params, :mu, nothing) : val
    end
    if mu === nothing
        mu_mev = _get_param(params, :mu_mev, nothing)
        mu = mu_mev === nothing ? nothing : mu_mev / ħc_MeV_fm
    end

    rho = _get_param(params, :rho, nothing)
    xi = _get_param(params, :xi, nothing)

    return (T = T, mu = mu, rho = rho, xi = xi)
end

function _distance(seed::SeedPoint, target, weights)
    total = 0.0
    used = 0
    for axis in (:T, :mu, :rho, :xi)
        val = getfield(target, axis)
        if val !== nothing
            weight = getfield(weights, axis)
            diff = (val - getfield(seed, axis)) * weight
            total += diff^2
            used += 1
        end
    end
    return used == 0 ? Inf : sqrt(total)
end

# ---------------- public API ----------------

function load_seed_table(; path::AbstractString = DEFAULT_SEED_PATH, force::Bool = false)
    if force || _seed_cache[] === nothing || _seed_cache_path[] != path
        if !isfile(path)
            error("Seed table not found: $path")
        end
        data, header = readdlm(path, ',', header = true)
        header_syms = Symbol.(vec(header))
        col_index = Dict(header_syms[i] => i for i in eachindex(header_syms))
        function col(row, name)
            row[col_index[name]]
        end
        seeds = SeedPoint[]
        for r in 1:size(data, 1)
            row = data[r, :]
            push!(seeds, SeedPoint(
                col(row, :T_MeV) / ħc_MeV_fm,
                col(row, :mu_MeV) / ħc_MeV_fm,
                col(row, :rho),
                col(row, :xi),
                col(row, :phi_u),
                col(row, :phi_d),
                col(row, :phi_s),
                col(row, :Phi1),
                col(row, :Phi2),
                col(row, :mu_u_MeV) / ħc_MeV_fm,
                col(row, :mu_d_MeV) / ħc_MeV_fm,
                col(row, :mu_s_MeV) / ħc_MeV_fm,
                path,
            ))
        end
        _seed_cache[] = seeds
        _seed_cache_path[] = path
    end
    return _seed_cache[]
end

function seed_state(seed::SeedPoint)
    return (
        x = [
            seed.phi_u,
            seed.phi_d,
            seed.phi_s,
            seed.Phi1,
            seed.Phi2,
            seed.mu_u,
            seed.mu_d,
            seed.mu_s,
        ],
        mu = (seed.mu_u, seed.mu_d, seed.mu_s),
        phi = (seed.phi_u, seed.phi_d, seed.phi_s),
        Polyakov = (seed.Phi1, seed.Phi2),
    )
end

function find_initial_seed(params; weights = DEFAULT_WEIGHTS, path::AbstractString = DEFAULT_SEED_PATH)
    target = _normalize_request(params)
    seeds = load_seed_table(path = path)
    best_seed = nothing
    best_dist = Inf
    for seed in seeds
        dist = _distance(seed, target, weights)
        if dist < best_dist
            best_seed = seed
            best_dist = dist
        end
    end
    if best_seed === nothing
        error("没有可用的种子点，请先提供至少一个种子参数")
    end
    return (
        seed = best_seed,
        distance = best_dist,
        used_axes = count(x -> x !== nothing, (target.T, target.mu, target.rho, target.xi)),
        state = seed_state(best_seed),
    )
end

end # module
