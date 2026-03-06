#!/usr/bin/env julia
using BenchmarkTools
using Random

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
Models.pnjl_module()

using .Constants_PNJL: ħc_MeV_fm
const PNJL = Models.pnjl_module()

const DEFAULT_STRATEGY = PNJL.DefaultSeed()
const MULTI_STRATEGY = PNJL.MultiSeed()

function random_request()
    return (
        T_mev = 40.0 + rand() * 180.0,
        mu_mev = -50.0 + rand() * 550.0,
        xi = rand(),
    )
end

function benchmark_once()
    request = random_request()
    T_fm = request.T_mev / ħc_MeV_fm
    mu_fm = request.mu_mev / ħc_MeV_fm
    theta = [T_fm, mu_fm]
    mode = PNJL.FixedMu()

    println("Benchmark target: T=$(request.T_mev) MeV, mu=$(request.mu_mev) MeV, xi=$(request.xi)")

    println("\n@btime DefaultSeed get_seed")
    @btime PNJL.get_seed($DEFAULT_STRATEGY, $theta, $mode)

    println("\n@btime MultiSeed get_seed (first candidate)")
    @btime PNJL.get_seed($MULTI_STRATEGY, $theta, $mode)

    println("\n@btime MultiSeed get_all_seeds")
    @btime PNJL.get_all_seeds($MULTI_STRATEGY, $theta, $mode)
end

benchmark_once()
