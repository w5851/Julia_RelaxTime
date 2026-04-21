#!/usr/bin/env julia

# Benchmark: PNJLModel.number_densities baseline
#
# Goals:
# - Measure the current cost/allocation profile of Models.number_densities.
# - Compare default cached-node retrieval with explicit thermal_nodes reuse.
# - Provide evidence on whether models-side node handling is still a material hotspot.
#
# Run:
#   julia --project=. benchmark/pnjl/bench_number_densities.jl

using BenchmarkTools
using Printf
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const MODEL = Models.PNJLModel()
const STATE = Models.MeanFieldState(SVector{3, Float64}(-1.5, -1.5, -2.0); Phi=0.3, PhiBar=0.3)
const TEMPERATURE = 0.5
const MU_VEC = SVector{3, Float64}(1.0, 1.0, 1.0)
const P_NUM = 24
const T_NUM = 6
const XI = 0.0
const THERMAL_NODES = Models.cached_nodes(P_NUM, T_NUM)

function summarize_trial(name::String, trial::BenchmarkTools.Trial)
    estimate = median(trial)
    return (
        name=name,
        time_ms=estimate.time / 1.0e6,
        memory_kib=estimate.memory / 1024.0,
        allocs=estimate.allocs,
    )
end

function print_summary(summary)
    @printf("%-44s %10.3f ms %10.1f KiB %8d allocs\n", summary.name, summary.time_ms, summary.memory_kib, summary.allocs)
end

function register_suite!()
    isdefined(Main, :SUITE) || return nothing

    Main.SUITE["pnjl"]["number_densities_default"] = @benchmarkable Models.number_densities(
        $MODEL,
        $STATE,
        $TEMPERATURE,
        $MU_VEC;
        p_num=$P_NUM,
        t_num=$T_NUM,
        xi=$XI,
    ) evals=1

    Main.SUITE["pnjl"]["number_densities_reuse_thermal_nodes"] = @benchmarkable Models.number_densities(
        $MODEL,
        $STATE,
        $TEMPERATURE,
        $MU_VEC;
        thermal_nodes=$THERMAL_NODES,
        xi=$XI,
    ) evals=1

    return nothing
end

function run_report(; samples::Int=20)
    bench_default = @benchmark Models.number_densities(
        $MODEL,
        $STATE,
        $TEMPERATURE,
        $MU_VEC;
        p_num=$P_NUM,
        t_num=$T_NUM,
        xi=$XI,
    ) samples=samples evals=1

    bench_reuse = @benchmark Models.number_densities(
        $MODEL,
        $STATE,
        $TEMPERATURE,
        $MU_VEC;
        thermal_nodes=$THERMAL_NODES,
        xi=$XI,
    ) samples=samples evals=1

    println("PNJLModel.number_densities baseline")
    println(repeat("=", 84))
    @printf("nodes=(p=%d, t=%d), xi=%.2f\n\n", P_NUM, T_NUM, XI)
    print_summary(summarize_trial("number_densities default", bench_default))
    print_summary(summarize_trial("number_densities reuse thermal_nodes", bench_reuse))

    ratio = median(bench_reuse).time / median(bench_default).time
    println("\nDecision hints")
    println(repeat("-", 84))
    @printf("reuse_thermal_nodes / default median time: %.2f%%\n", 100 * ratio)
    return nothing
end

register_suite!()

if abspath(PROGRAM_FILE) == @__FILE__
    samples = parse(Int, get(ENV, "BENCH_SAMPLES_PNJL_DENSITY", "20"))
    run_report(; samples=samples)
end
