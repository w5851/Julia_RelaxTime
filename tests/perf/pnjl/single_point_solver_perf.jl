# Benchmark a single PNJL gap solve under varying NLsolve configurations (using new architecture)
using Printf
using Statistics
using BenchmarkTools
using JSON3
using Dates
using NLsolve: LineSearches

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PERF_OUTPUT_DIR = joinpath(PROJECT_ROOT, "tests", "perf", "results", "pnjl")
const JSON_OUTPUT = joinpath(PERF_OUTPUT_DIR, "single_point_benchmark.json")
const MARKDOWN_OUTPUT = joinpath(PERF_OUTPUT_DIR, "single_point_benchmark.md")

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .PNJL: solve, FixedRho

const BENCH_CONFIG = (
    T_mev = 50.0,
    rho_target = 0.05,
    xi = 0.0,
    p_num = 24,
    t_num = 12,
)
struct BenchResult
    label::String
    params::String
    stats::NamedTuple
    trial::BenchmarkTools.Trial
end


const DEFAULT_SAMPLES = 20

struct SolverBenchmark
    label::String
    kwargs::NamedTuple
end

const BENCHMARKS = [
    SolverBenchmark("Newton (default line search)", (method = :newton,)),
    SolverBenchmark(
        "Newton + BackTracking line search",
        (method = :newton, linesearch = LineSearches.BackTracking()),
    ),
    SolverBenchmark("Trust-region", (method = :trust_region,)),
]

function run_single_point(; kwargs...)
    T_fm = BENCH_CONFIG.T_mev / ħc_MeV_fm
    return solve(FixedRho(BENCH_CONFIG.rho_target), T_fm;
        xi = BENCH_CONFIG.xi,
        p_num = BENCH_CONFIG.p_num,
        t_num = BENCH_CONFIG.t_num,
        kwargs...,
    )
end

run_single_point_named(named_kwargs::NamedTuple) = run_single_point(; named_kwargs...)

function benchmark_solver(named_kwargs; samples::Int = DEFAULT_SAMPLES)
    run_single_point_named(named_kwargs) # warm-up JIT/cache
    local_kwargs = named_kwargs
    bench = @benchmarkable run_single_point_named($local_kwargs) evals=1
    return run(bench; samples=samples)
end

function summarize_trial(trial::BenchmarkTools.Trial)
    times_ms = trial.times ./ 1.0e6
    return (;
        samples = length(times_ms),
        min = minimum(times_ms),
        median = median(times_ms),
        max = maximum(times_ms),
        mean = mean(times_ms),
        allocs = trial.allocs,
        memory_bytes = trial.memory,
        gctime_ms = mean(trial.gctimes) / 1.0e6,
    )
end

function write_reports(results::Vector{BenchResult})
    mkpath(PERF_OUTPUT_DIR)
    timestamp = Dates.format(Dates.now(), Dates.ISODateTimeFormat)
    json_payload = (
        generated_at = timestamp,
        samples = DEFAULT_SAMPLES,
        config = BENCH_CONFIG,
        benchmarks = [
            (
                label = r.label,
                params = r.params,
                stats = r.stats,
            ) for r in results
        ],
    )
    open(JSON_OUTPUT, "w") do io
        JSON3.pretty(io, json_payload)
    end

    open(MARKDOWN_OUTPUT, "w") do io
        println(io, "# PNJL Single-Point Benchmark")
        println(io, "Generated at: $timestamp")
        println(io, "")
        println(io, "| label | params | min (ms) | median (ms) | mean (ms) | max (ms) |")
        println(io, "| --- | --- | ---: | ---: | ---: | ---: |")
        for r in results
            s = r.stats
            println(io, "| $(r.label) | $(r.params) | $(round(s.min, digits=3)) | $(round(s.median, digits=3)) | $(round(s.mean, digits=3)) | $(round(s.max, digits=3)) |")
        end
    end

    println("Saved benchmark reports to:\n  JSON -> $JSON_OUTPUT\n  Markdown -> $MARKDOWN_OUTPUT")
end

function main()
    @printf("PNJL single-point benchmark: T=%.2f MeV, rho=%.4f, xi=%.2f (p=%d, t=%d)\n",
        BENCH_CONFIG.T_mev,
        BENCH_CONFIG.rho_target,
        BENCH_CONFIG.xi,
        BENCH_CONFIG.p_num,
        BENCH_CONFIG.t_num,
    )
    results = BenchResult[]
    @printf("BenchmarkTools samples per configuration: %d (evals=1)\n\n", DEFAULT_SAMPLES)
    for bench in BENCHMARKS
        println("=== ", bench.label, " ===")
        trial = benchmark_solver(bench.kwargs)
        stats = summarize_trial(trial)
        @printf("  min   : %8.3f ms\n", stats.min)
        @printf("  median: %8.3f ms\n", stats.median)
        @printf("  mean  : %8.3f ms\n", stats.mean)
        @printf("  max   : %8.3f ms\n\n", stats.max)
        push!(results, BenchResult(bench.label, string(bench.kwargs), stats, trial))
    end
    write_reports(results)
end

main()
