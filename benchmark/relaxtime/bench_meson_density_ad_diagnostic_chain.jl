#!/usr/bin/env julia

# Benchmark: meson-density AD diagnostic chain
#
# Goals:
# - Measure absolute cost of the workflow-level AD diagnostic point.
# - Measure relative cost share against equilibrium solve and meson-density point workflow.
# - Provide a stable baseline before promoting AD as production dδ/dω derivative path.
#
# Run:
#   julia --project=benchmark benchmark/relaxtime/bench_meson_density_ad_diagnostic_chain.jl
#
# Optional env knobs:
#   BENCH_MD_AD_T_MEV=210.0
#   BENCH_MD_AD_XI=0.0
#   BENCH_MD_AD_Q=0.5
#   BENCH_MD_AD_OMEGA=0.6
#   BENCH_MD_AD_FD_STEP=1e-5
#   BENCH_MD_AD_EQ_SAMPLES=8
#   BENCH_MD_AD_DIAG_SAMPLES=8
#   BENCH_MD_AD_DENSITY_SAMPLES=6
#   BENCH_MD_AD_FULL_SAMPLES=5

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const BENCH_ENV = normpath(joinpath(@__DIR__, ".."))

pushfirst!(LOAD_PATH, PROJECT_ROOT)
pushfirst!(LOAD_PATH, BENCH_ENV)

using BenchmarkTools
using Printf

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .Models: solve_gap_and_meson_point,
               solve_phase_shift_point_diagnostic_from_meson_point,
               solve_phase_shift_meson_density_from_meson_point,
               solve_gap_and_phase_shift_meson_density_point

const T_MEV = parse(Float64, get(ENV, "BENCH_MD_AD_T_MEV", "210.0"))
const XI = parse(Float64, get(ENV, "BENCH_MD_AD_XI", "0.0"))
const Q0 = parse(Float64, get(ENV, "BENCH_MD_AD_Q", "0.5"))
const OMEGA0 = parse(Float64, get(ENV, "BENCH_MD_AD_OMEGA", "0.6"))
const FD_STEP = parse(Float64, get(ENV, "BENCH_MD_AD_FD_STEP", "1e-5"))

const EQ_SAMPLES = parse(Int, get(ENV, "BENCH_MD_AD_EQ_SAMPLES", "8"))
const DIAG_SAMPLES = parse(Int, get(ENV, "BENCH_MD_AD_DIAG_SAMPLES", "8"))
const DENSITY_SAMPLES = parse(Int, get(ENV, "BENCH_MD_AD_DENSITY_SAMPLES", "6"))
const FULL_SAMPLES = parse(Int, get(ENV, "BENCH_MD_AD_FULL_SAMPLES", "5"))

const T_FM = T_MEV / ħc_MeV_fm

function build_equilibrium_point()
    return solve_gap_and_meson_point(
        T_FM,
        0.0;
        xi=XI,
        p_num=12,
        t_num=6,
    )
end

const EQ_POINT = build_equilibrium_point()

function ad_diag_post_equilibrium()
    return solve_phase_shift_point_diagnostic_from_meson_point(
        EQ_POINT;
        mesons=(:pi,),
        q_values=[Q0],
        omega_values=[OMEGA0],
        scheme=:current,
        fd_step=FD_STEP,
    )
end

function density_post_equilibrium()
    return solve_phase_shift_meson_density_from_meson_point(
        EQ_POINT;
        scheme=:current,
    )
end

function density_full_point()
    return solve_gap_and_phase_shift_meson_density_point(
        T_FM,
        0.0;
        xi=XI,
        p_num=12,
        t_num=6,
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
    @printf("%-38s %10.3f ms %10.1f KiB %8d allocs\n",
        summary.name, summary.time_ms, summary.memory_kib, summary.allocs)
end

function register_suite!()
    isdefined(Main, :SUITE) || return nothing
    Main.SUITE["relaxtime"]["meson_density_equilibrium_point"] =
        @benchmarkable build_equilibrium_point() evals=1
    Main.SUITE["relaxtime"]["meson_density_ad_diag_post_eq"] =
        @benchmarkable ad_diag_post_equilibrium() evals=1
    Main.SUITE["relaxtime"]["meson_density_phase_shift_post_eq"] =
        @benchmarkable density_post_equilibrium() evals=1
    Main.SUITE["relaxtime"]["meson_density_phase_shift_full_point"] =
        @benchmarkable density_full_point() evals=1
    return nothing
end

function run_report()
    # warmup
    ad_diag_post_equilibrium()
    density_post_equilibrium()
    density_full_point()

    eq_trial = @benchmark build_equilibrium_point() samples=EQ_SAMPLES evals=1
    diag_trial = @benchmark ad_diag_post_equilibrium() samples=DIAG_SAMPLES evals=1
    density_trial = @benchmark density_post_equilibrium() samples=DENSITY_SAMPLES evals=1
    full_trial = @benchmark density_full_point() samples=FULL_SAMPLES evals=1

    println("Meson-density AD diagnostic benchmark")
    println(repeat("=", 88))
    @printf("T=%.1f MeV, xi=%.3f, q=%.3f, omega=%.3f, fd_step=%g\n\n",
        T_MEV, XI, Q0, OMEGA0, FD_STEP)

    print_summary(summarize_trial("equilibrium meson point", eq_trial))
    print_summary(summarize_trial("AD diagnostic post-equilibrium", diag_trial))
    print_summary(summarize_trial("phase-shift density post-equilibrium", density_trial))
    print_summary(summarize_trial("full phase-shift density point", full_trial))

    eq_t = median(eq_trial).time
    diag_t = median(diag_trial).time
    density_t = median(density_trial).time
    full_t = median(full_trial).time

    println("\nRelative share")
    println(repeat("-", 88))
    @printf("AD diagnostic / equilibrium point: %.2f%%\n", 100 * diag_t / eq_t)
    @printf("AD diagnostic / full density point: %.2f%%\n", 100 * diag_t / full_t)
    @printf("post-eq density / full density point: %.2f%%\n", 100 * density_t / full_t)
    @printf("equilibrium / full density point: %.2f%%\n", 100 * eq_t / full_t)
    return nothing
end

register_suite!()

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_report()
end
