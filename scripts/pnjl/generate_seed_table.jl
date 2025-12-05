#!/usr/bin/env julia
# Generate PNJL seed table via Sobol/LHS sampling for warm-start solvers

using Random
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))

module PNJLLocal
    using ..Constants_PNJL
    include(joinpath(Main.PROJECT_ROOT, "src", "pnjl", "solvers", "AnisoGapSolver.jl"))
    include(joinpath(Main.PROJECT_ROOT, "src", "pnjl", "seeds", "TrhoSeedChain.jl"))
    export AnisoGapSolver, TrhoSeedChain
end

using .PNJLLocal: AnisoGapSolver, TrhoSeedChain
const DEFAULT_MU_GUESS = AnisoGapSolver.DEFAULT_MU_GUESS
const ħc_MeV_fm = Main.Constants_PNJL.ħc_MeV_fm
const SEED_T_DIGITS = 1
const SEED_XI_DIGITS = 2
const ACCEPTABLE_RESIDUAL = 1e-4
const DEFAULT_CHAIN_DIR = TrhoSeedChain.DEFAULT_CHAIN_DIR

const OUTPUT_DIR = joinpath(PROJECT_ROOT, "data", "raw", "pnjl", "seeds")
const DEFAULT_OUTPUT = joinpath(OUTPUT_DIR, "sobol_seed_table.csv")
const HEADER = [
    "T_MeV",
    "mu_MeV",
    "rho",
    "xi",
    "phi_u",
    "phi_d",
    "phi_s",
    "Phi1",
    "Phi2",
    "mu_u_MeV",
    "mu_d_MeV",
    "mu_s_MeV",
]

const DEFAULT_CONFIG = (
    samples = 160,
    T_range = (1.0, 400.0),
    mu_range = (0.0, 500.0),
    xi_range = (-0.8, 0.8),
    rng_seed = 42,
    p_num = 24,
    t_num = 8,
)

function parse_args()
    cfg = Dict{Symbol, Any}(
        :samples => DEFAULT_CONFIG.samples,
        :output => DEFAULT_OUTPUT,
        :rng_seed => DEFAULT_CONFIG.rng_seed,
        :p_num => DEFAULT_CONFIG.p_num,
        :t_num => DEFAULT_CONFIG.t_num,
        :T_range => DEFAULT_CONFIG.T_range,
        :mu_range => DEFAULT_CONFIG.mu_range,
        :xi_range => DEFAULT_CONFIG.xi_range,
        :chain_dir => DEFAULT_CHAIN_DIR,
        :use_chain => true,
        :chain_t_tol => 3.0,
        :chain_xi_tol => 1e-3,
    )
    for arg in ARGS
        if startswith(arg, "--")
            key, value = split(arg[3:end], "=", limit = 2)
            sym = Symbol(key)
            if sym == :output
                cfg[:output] = value
            elseif sym == :samples
                cfg[:samples] = parse(Int, value)
            elseif sym == :rng_seed
                cfg[:rng_seed] = parse(Int, value)
            elseif sym == :p_num
                cfg[:p_num] = parse(Int, value)
            elseif sym == :t_num
                cfg[:t_num] = parse(Int, value)
            elseif sym == :T_min
                lo, hi = cfg[:T_range]
                cfg[:T_range] = (parse(Float64, value), hi)
            elseif sym == :T_max
                lo, _ = cfg[:T_range]
                cfg[:T_range] = (lo, parse(Float64, value))
            elseif sym == :mu_min
                lo, hi = cfg[:mu_range]
                cfg[:mu_range] = (parse(Float64, value), hi)
            elseif sym == :mu_max
                lo, _ = cfg[:mu_range]
                cfg[:mu_range] = (lo, parse(Float64, value))
            elseif sym == :xi_max
                cfg[:xi_range] = (cfg[:xi_range][1], parse(Float64, value))
            elseif sym == :xi_min
                cfg[:xi_range] = (parse(Float64, value), cfg[:xi_range][2])
            elseif sym == :chain_dir
                cfg[:chain_dir] = value
            elseif sym == :disable_chain
                cfg[:use_chain] = false
            elseif sym == :chain_t_tol
                cfg[:chain_t_tol] = parse(Float64, value)
            elseif sym == :chain_xi_tol
                cfg[:chain_xi_tol] = parse(Float64, value)
            else
                @warn "Unknown argument" arg
            end
        end
    end
    return cfg
end

function latin_hypercube(n::Int, d::Int; rng = Random.default_rng())
    samples = zeros(n, d)
    for dim in 1:d
        perm = randperm(rng, n)
        for i in 1:n
            samples[i, dim] = (perm[i] - rand(rng)) / n
        end
    end
    return samples
end

scale(val, range::Tuple{Float64, Float64}) = range[1] + val * (range[2] - range[1])

function sample_parameters(cfg)
    rng = MersenneTwister(cfg[:rng_seed])
    base = latin_hypercube(cfg[:samples], 3; rng = rng)
    params = Vector{NamedTuple{(:T_mev, :mu_mev, :xi), NTuple{3, Float64}}}()
    for i in 1:cfg[:samples]
        push!(params, (
            T_mev = scale(base[i, 1], cfg[:T_range]),
            mu_mev = scale(base[i, 2], cfg[:mu_range]),
            xi = scale(base[i, 3], cfg[:xi_range]),
        ))
    end
    sort!(params, by = p -> (_seed_key(p.T_mev, p.xi), p.mu_mev))
    return params
end

function build_row(T_mev, mu_mev, xi, result)
    phi_u, phi_d, phi_s = result.solution[1:3]
    Phi1, Phi2 = result.solution[4], result.solution[5]
    mu_u = result.mu_vec[1] * ħc_MeV_fm
    mu_d = result.mu_vec[2] * ħc_MeV_fm
    mu_s = result.mu_vec[3] * ħc_MeV_fm
    return (
        T_mev,
        mu_mev,
        result.rho,
        xi,
        phi_u,
        phi_d,
        phi_s,
        Phi1,
        Phi2,
        mu_u,
        mu_d,
        mu_s,
    )
end

function format_row(row)
    return join(string.(row), ",")
end

function main()
    cfg = parse_args()
    params = sample_parameters(cfg)
    mkpath(dirname(cfg[:output]))
    successes = 0
    failures = 0
    continuation = Dict{Tuple{Float64, Float64}, Vector{Float64}}()
    io = open(cfg[:output], "w")
    println(io, join(HEADER, ","))
    for (idx, p) in enumerate(params)
        mu_fm = p.mu_mev / ħc_MeV_fm
        seed, _ = _select_seed(p, continuation, cfg)
        try
            res = AnisoGapSolver.solve_fixed_mu(p.T_mev, mu_fm;
                xi = p.xi,
                p_num = cfg[:p_num],
                t_num = cfg[:t_num],
                seed_state = seed,
            )
            res, refined = _refine_if_needed(p.T_mev, mu_fm, p.xi, res; p_num=cfg[:p_num], t_num=cfg[:t_num])
            if res !== nothing && res.converged
                row = format_row(build_row(p.T_mev, p.mu_mev, p.xi, res))
                println(io, row)
                continuation[_seed_key(p.T_mev, p.xi)] = copy(res.solution)
                successes += 1
            else
                failures += 1
            end
        catch err
            @warn "Seed solve failed" idx err
            failures += 1
        end
    end
    close(io)
    @printf("Generated %d seeds (failures=%d) -> %s\n", successes, failures, cfg[:output])
end

_seed_key(T, xi) = (round(Float64(T); digits=SEED_T_DIGITS), round(Float64(xi); digits=SEED_XI_DIGITS))

function _select_seed(params, continuation, cfg)
    if cfg[:use_chain]
        seed = TrhoSeedChain.find_chain_seed(params.T_mev, params.xi, params.mu_mev;
            chain_dir=cfg[:chain_dir], T_tol=cfg[:chain_t_tol], xi_tol=cfg[:chain_xi_tol])
        if seed !== nothing
            return seed, "trho-chain"
        end
    end
    key = _seed_key(params.T_mev, params.xi)
    if haskey(continuation, key)
        return copy(continuation[key]), "continuation"
    end
    return copy(DEFAULT_MU_GUESS), "default"
end

function _refine_if_needed(T_mev, mu_fm, xi, result; p_num, t_num)
    result === nothing && return nothing, false
    if result.converged
        return result, false
    end
    residual = result.residual_norm
    if !isfinite(residual) || residual > ACCEPTABLE_RESIDUAL
        return result, false
    end
    try
        refined = AnisoGapSolver.solve_fixed_mu(T_mev, mu_fm;
            xi = xi,
            p_num = p_num,
            t_num = t_num,
            seed_state = result.solution,
        )
        return refined, true
    catch
        return result, false
    end
end

main()
