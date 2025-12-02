#!/usr/bin/env julia
using BenchmarkTools
using Random

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "SeedCache.jl"))

const WEIGHTS = SeedCache.DEFAULT_WEIGHTS
const PATH = SeedCache.DEFAULT_SEED_PATH
const SEEDS = SeedCache.load_seed_table(path = PATH)

function random_request()
    return Dict(
        :T_mev => 40.0 + rand() * 180.0,
        :mu_mev => -50.0 + rand() * 550.0,
        :xi => rand(),
    )
end

function benchmark_once()
    request = random_request()
    target = SeedCache._normalize_request(request)
    neighbors = SeedCache._query_mu_neighbors(SEEDS, target, WEIGHTS, PATH, 3)
    println("Benchmark target: T=$(target.T*Constants_PNJL.ħc_MeV_fm) MeV, mu=$(target.mu*Constants_PNJL.ħc_MeV_fm) MeV, xi=$(target.xi)")
    println("Near neighbor distances: ", [n.distance for n in neighbors])

    println("\n@btime find_initial_seed (k=3)")
    @btime SeedCache.find_initial_seed($request; weights = WEIGHTS, path = PATH, k_neighbors = 3)

    println("\n@btime neighbor query only")
    @btime SeedCache._query_mu_neighbors($SEEDS, $target, WEIGHTS, PATH, 3)

    println("\n@btime blend neighbors")
    @btime SeedCache._blend_neighbors($neighbors)

    println("\n@btime find_initial_seed (k=1 no blend)")
    @btime SeedCache.find_initial_seed($request; weights = WEIGHTS, path = PATH, k_neighbors = 1)
end

benchmark_once()
