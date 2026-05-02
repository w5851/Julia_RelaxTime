#!/usr/bin/env julia

# Benchmark: Phase E3 BU double-integral diagnostic chain
#
# Goals:
# - Measure absolute cost of the current minimal BU double-integral point evaluation.
# - Separate equilibrium/meson-point solve from the post-equilibrium double integral.
# - Report relative cost share so later workflow/module promotion has a baseline.
#
# Run:
#   julia --project=benchmark benchmark/relaxtime/bench_bu_double_integral_phase_e3.jl
#
# Optional env knobs:
#   BENCH_BU_T_MEV=210.0
#   BENCH_BU_Q_MAX=12.0
#   BENCH_BU_Q_NODES=48
#   BENCH_BU_OMEGA_MAX=10.0
#   BENCH_BU_OMEGA_NODES=48
#   BENCH_BU_SAMPLES_EQ=8
#   BENCH_BU_SAMPLES_POST=8
#   BENCH_BU_SAMPLES_FULL=5

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const BENCH_ENV = normpath(joinpath(@__DIR__, ".."))

pushfirst!(LOAD_PATH, PROJECT_ROOT)
pushfirst!(LOAD_PATH, BENCH_ENV)

using BenchmarkTools
using Printf

const SCRIPT_PATH = joinpath(
    PROJECT_ROOT, "scripts", "analysis", "relaxtime", "diagnose_bu_double_integral_minimal.jl")
Base.include(Main, SCRIPT_PATH)

const T_MeV = parse(Float64, get(ENV, "BENCH_BU_T_MEV", "210.0"))
const Q_MAX = parse(Float64, get(ENV, "BENCH_BU_Q_MAX", "12.0"))
const Q_NODES = parse(Int, get(ENV, "BENCH_BU_Q_NODES", "48"))
const OMEGA_MAX = parse(Float64, get(ENV, "BENCH_BU_OMEGA_MAX", "10.0"))
const OMEGA_NODES = parse(Int, get(ENV, "BENCH_BU_OMEGA_NODES", "48"))
const ETA = parse(Float64, get(ENV, "BENCH_BU_ETA", "1e-6"))

const EQ_SAMPLES = parse(Int, get(ENV, "BENCH_BU_SAMPLES_EQ", "8"))
const POST_SAMPLES = parse(Int, get(ENV, "BENCH_BU_SAMPLES_POST", "8"))
const FULL_SAMPLES = parse(Int, get(ENV, "BENCH_BU_SAMPLES_FULL", "5"))

const Q_GRID_Q_W = Main.GaussLegendre.gauleg(0.0, Q_MAX, Q_NODES)
const O_GRID_O_W = Main.GaussLegendre.gauleg(0.05, OMEGA_MAX, OMEGA_NODES)
const Q_GRID = Q_GRID_Q_W[1]
const Q_W = Q_GRID_Q_W[2]
const O_GRID = O_GRID_O_W[1]
const O_W = O_GRID_O_W[2]
const EQ_POINT = Main.solve_equilibrium_meson_point(T_MeV)

function bu_full_point()
    res = Main.solve_equilibrium_meson_point(T_MeV)
    return Main.compute_double_integral_rows(
        res;
        q_nodes=Q_GRID,
        q_weights=Q_W,
        omega_nodes=O_GRID,
        omega_weights=O_W,
        qmax_request=Q_MAX,
        omega_max_request=OMEGA_MAX,
        eta=ETA,
    )
end

function bu_post_equilibrium()
    return Main.compute_double_integral_rows(
        EQ_POINT;
        q_nodes=Q_GRID,
        q_weights=Q_W,
        omega_nodes=O_GRID,
        omega_weights=O_W,
        qmax_request=Q_MAX,
        omega_max_request=OMEGA_MAX,
        eta=ETA,
    )
end

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
    @printf("%-42s %10.3f ms %10.1f KiB %8d allocs\n",
        summary.name, summary.time_ms, summary.memory_kib, summary.allocs)
end

function register_suite!()
    isdefined(Main, :SUITE) || return nothing
    Main.SUITE["relaxtime"]["bu_phase_e3_equilibrium_point"] =
        @benchmarkable Main.solve_equilibrium_meson_point($T_MeV) evals=1
    Main.SUITE["relaxtime"]["bu_phase_e3_post_equilibrium"] =
        @benchmarkable bu_post_equilibrium() evals=1
    Main.SUITE["relaxtime"]["bu_phase_e3_full_point"] =
        @benchmarkable bu_full_point() evals=1
    return nothing
end

function run_report()
    eq_trial = @benchmark Main.solve_equilibrium_meson_point($T_MeV) samples=EQ_SAMPLES evals=1
    post_trial = @benchmark bu_post_equilibrium() samples=POST_SAMPLES evals=1
    full_trial = @benchmark bu_full_point() samples=FULL_SAMPLES evals=1

    println("Phase E3 BU double-integral benchmark")
    println(repeat("=", 88))
    @printf("T=%.1f MeV, q_max=%.1f, q_nodes=%d, omega_max=%.1f, omega_nodes=%d\n\n",
        T_MeV, Q_MAX, Q_NODES, OMEGA_MAX, OMEGA_NODES)

    print_summary(summarize_trial("equilibrium meson point", eq_trial))
    print_summary(summarize_trial("post-equilibrium double integral", post_trial))
    print_summary(summarize_trial("full BU point", full_trial))

    eq_t = median(eq_trial).time
    post_t = median(post_trial).time
    full_t = median(full_trial).time

    println("\nRelative share")
    println(repeat("-", 88))
    @printf("equilibrium / full point: %.2f%%\n", 100 * eq_t / full_t)
    @printf("post-equilibrium integral / full point: %.2f%%\n", 100 * post_t / full_t)
    return nothing
end

register_suite!()

if abspath(PROGRAM_FILE) == @__FILE__
    run_report()
end
