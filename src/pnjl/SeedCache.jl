module SeedCache

using DelimitedFiles
using NearestNeighbors
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

const _SEED_DIR = normpath(joinpath(@__DIR__, "..", "..", "data", "raw", "pnjl", "seeds"))
const GENERATED_SEED_FILE = "sobol_seed_table.csv"
const DEFAULT_SEED_PATH = normpath(joinpath(_SEED_DIR, GENERATED_SEED_FILE))
const FALLBACK_SEED_PATH = normpath(joinpath(_SEED_DIR, "aniso_coarse.csv"))
const DEFAULT_WEIGHTS = (T = 1.0, mu = 1.0, rho = 0.6, xi = 0.4)

const _seed_cache = Ref(Vector{SeedPoint}())
const _seed_cache_path = Ref{String}("")
const _mu_tree_cache = Ref{Union{Nothing, Tuple{KDTree, String, NTuple{3, Float64}}}}(nothing)
const _rho_tree_cache = Ref{Union{Nothing, Tuple{KDTree, String, NTuple{3, Float64}}}}(nothing)

_weights_key_mu(w) = (w.T, w.mu, w.xi)
_weights_key_rho(w) = (w.T, w.rho, w.xi)

function _resolve_seed_path(path::AbstractString)
    if isfile(path)
        return path
    elseif path == DEFAULT_SEED_PATH && isfile(FALLBACK_SEED_PATH)
        @warn "Seed table $path not found, falling back to $FALLBACK_SEED_PATH"
        return FALLBACK_SEED_PATH
    else
        error("Seed table not found: $path")
    end
end

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

function _collect_neighbors(seeds, indices, distances)
    count = length(indices)
    neighbors = Vector{NamedTuple}(undef, count)
    for j in 1:count
        seed = seeds[indices[j]]
        state = seed_state(seed)
        neighbors[j] = (seed = seed, distance = distances[j], state = state)
    end
    return sort(neighbors, by = n -> n.distance)
end

function _linear_neighbors(seeds, target, weights, k)
    scored = NamedTuple[]
    for seed in seeds
        dist = _distance(seed, target, weights)
        push!(scored, (seed = seed, distance = dist, state = seed_state(seed)))
    end
    sort!(scored, by = n -> n.distance)
    return scored[1:clamp(k, 1, length(scored))]
end

function _blend_neighbors(neighbors)
    isempty(neighbors) && return nothing
    weights = similar(neighbors, Float64)
    for (i, neighbor) in enumerate(neighbors)
        weights[i] = neighbor.distance == 0.0 ? Inf : 1.0 / neighbor.distance
    end
    if any(isinf, weights)
        idx = findfirst(isinf, weights)
        rho_seed = copy(neighbors[idx].state.x)
        mu_seed = rho_seed[1:5]
        return (rho_seed = rho_seed, mu_seed = mu_seed)
    end
    norm = sum(weights)
    rho_seed = zeros(Float64, length(neighbors[1].state.x))
    mu_seed = zeros(Float64, 5)
    for (neighbor, w) in zip(neighbors, weights)
        rho_seed .+= w .* neighbor.state.x
        mu_seed .+= w .* neighbor.state.x[1:5]
    end
    rho_seed ./= norm
    mu_seed ./= norm
    return (rho_seed = rho_seed, mu_seed = mu_seed)
end

function _query_mu_neighbors(seeds, target, weights, path, k)
    tree = _ensure_mu_tree(weights, path)
    xi_val = target.xi === nothing ? 0.0 : target.xi
    point = [weights.T * target.T, weights.mu * target.mu, weights.xi * xi_val]
    idxs, dists = knn(tree, point, clamp(k, 1, length(seeds)))
    return _collect_neighbors(seeds, idxs, dists)
end

function _query_rho_neighbors(seeds, target, weights, path, k)
    tree = _ensure_rho_tree(weights, path)
    xi_val = target.xi === nothing ? 0.0 : target.xi
    point = [weights.T * target.T, weights.rho * target.rho, weights.xi * xi_val]
    idxs, dists = knn(tree, point, clamp(k, 1, length(seeds)))
    return _collect_neighbors(seeds, idxs, dists)
end

function _build_mu_tree(seeds, weights, path)
    coords = zeros(3, length(seeds))
    for (i, seed) in enumerate(seeds)
        coords[1, i] = weights.T * seed.T
        coords[2, i] = weights.mu * seed.mu
        coords[3, i] = weights.xi * seed.xi
    end
    tree = KDTree(coords)
    _mu_tree_cache[] = (tree, path, _weights_key_mu(weights))
    return tree
end

function _build_rho_tree(seeds, weights, path)
    coords = zeros(3, length(seeds))
    for (i, seed) in enumerate(seeds)
        coords[1, i] = weights.T * seed.T
        coords[2, i] = weights.rho * seed.rho
        coords[3, i] = weights.xi * seed.xi
    end
    tree = KDTree(coords)
    _rho_tree_cache[] = (tree, path, _weights_key_rho(weights))
    return tree
end

function _ensure_mu_tree(weights, path)
    resolved = _resolve_seed_path(path)
    cache = _mu_tree_cache[]
    if cache !== nothing && cache[2] == resolved && cache[3] == _weights_key_mu(weights)
        return cache[1]
    end
    seeds = load_seed_table(path = resolved)
    return _build_mu_tree(seeds, weights, resolved)
end

function _ensure_rho_tree(weights, path)
    resolved = _resolve_seed_path(path)
    cache = _rho_tree_cache[]
    if cache !== nothing && cache[2] == resolved && cache[3] == _weights_key_rho(weights)
        return cache[1]
    end
    seeds = load_seed_table(path = resolved)
    return _build_rho_tree(seeds, weights, resolved)
end

# ---------------- public API ----------------

function load_seed_table(; path::AbstractString = DEFAULT_SEED_PATH, force::Bool = false)
    resolved = _resolve_seed_path(path)
    if force || _seed_cache[] === nothing || _seed_cache_path[] != resolved
        data, header = readdlm(resolved, ',', header = true)
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
                resolved,
            ))
        end
        _seed_cache[] = seeds
        _seed_cache_path[] = resolved
        _mu_tree_cache[] = nothing
        _rho_tree_cache[] = nothing
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

function find_initial_seed(params; weights = DEFAULT_WEIGHTS, path::AbstractString = DEFAULT_SEED_PATH, k_neighbors::Int = 3)
    target = _normalize_request(params)
    target.T === nothing && error("seed lookup needs temperature T")
    seeds = load_seed_table(path = path)
    isempty(seeds) && error("seed table is empty")

    neighbors = if target.mu !== nothing
        _query_mu_neighbors(seeds, target, weights, path, k_neighbors)
    elseif target.rho !== nothing
        _query_rho_neighbors(seeds, target, weights, path, k_neighbors)
    else
        _linear_neighbors(seeds, target, weights, k_neighbors)
    end

    neighbors = isempty(neighbors) ? _linear_neighbors(seeds, target, weights, k_neighbors) : neighbors
    best = first(neighbors)
    blended = _blend_neighbors(neighbors)
    return (
        seed = best.seed,
        distance = best.distance,
        used_axes = count(x -> x !== nothing, (target.T, target.mu, target.rho, target.xi)),
        state = best.state,
        neighbors = neighbors,
        blended_state = blended,
    )
end

end # module
