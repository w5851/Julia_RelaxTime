#!/usr/bin/env julia

# Benchmark: number_density and average_scattering_rate baseline
#
# Goals:
# - Establish a comparable baseline for the density + average-scattering-rate chain.
# - Separate full-call cost from the semi-infinite momentum-grid mapping kernel.
# - Provide evidence for whether optimization should focus on cached grids or preallocation.
#
# Run:
#   julia --project=. benchmark/relaxtime/bench_density_and_scattering_rate.jl
#
# Optional env knobs:
#   BENCH_SAMPLES_DENSITY=20
#   BENCH_SAMPLES_ASR=5
#   BENCH_SAMPLES_KERNEL=200

using BenchmarkTools
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))

if !isdefined(Main, :RelaxTime)
    include(joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxTime.jl"))
end

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS, gauleg
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings

const ASR = AverageScatteringRate

function build_params(; T_MeV::Float64=150.0, muB_MeV::Float64=800.0,
    mu_s_MeV::Float64=0.0, m_u_MeV::Float64=300.0, m_d_MeV::Float64=300.0, m_s_MeV::Float64=500.0,
    phi::Float64=0.5, phibar::Float64=0.5, xi::Float64=0.8)

    T = T_MeV / ħc_MeV_fm
    μ_q = (muB_MeV / 3.0) / ħc_MeV_fm
    μ_u = μ_q
    μ_d = μ_q
    μ_s = mu_s_MeV / ħc_MeV_fm

    m_u = m_u_MeV / ħc_MeV_fm
    m_d = m_d_MeV / ħc_MeV_fm
    m_s = m_s_MeV / ħc_MeV_fm

    A_u = A(m_u, μ_u, T, phi, phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_d = A(m_d, μ_d, T, phi, phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_s = A(m_s, μ_s, T, phi, phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)

    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)

    quark_params = (m=(u=m_u, d=m_d, s=m_s), μ=(u=μ_u, d=μ_d, s=μ_s), A=(u=A_u, d=A_d, s=A_s))
    thermo_params = (T=T, Φ=phi, Φbar=phibar, ξ=xi)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    return (quark_params=quark_params, thermo_params=thermo_params, K_coeffs=K_coeffs)
end

function constant_sigma_cache(process::Symbol; sigma::Float64=1.0)
    cache = ASR.CrossSectionCache(process)
    ASR.insert_sigma!(cache, 0.0, sigma)
    ASR.insert_sigma!(cache, 500.0, sigma)
    return cache
end

function build_semi_infinite_momentum_grid(p_nodes::Int, scale::Float64)
    p_vals = Float64[]
    p_wts = Float64[]
    dp_jac = Float64[]
    t_grid, t_w = gauleg(0.0, 1.0, p_nodes)
    for (t, wt) in zip(t_grid, t_w)
        if t >= 0.9999
            continue
        end
        p = scale * t / (1.0 - t)
        dp_dt = scale / (1.0 - t)^2
        push!(p_vals, p)
        push!(p_wts, wt)
        push!(dp_jac, dp_dt)
    end
    return (p_vals=p_vals, p_wts=p_wts, dp_jac=dp_jac)
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
    @printf("%-40s %10.3f ms %10.1f KiB %8d allocs\n", summary.name, summary.time_ms, summary.memory_kib, summary.allocs)
end

const BENCH_PROCESS = :udbar_to_udbar
const BENCH_PARAMS = build_params()
const BENCH_CACHE = constant_sigma_cache(BENCH_PROCESS; sigma=1.0)
const DENSITY_FLAVOR = :u
const DENSITY_MASS = BENCH_PARAMS.quark_params.m.u
const DENSITY_MU = BENCH_PARAMS.quark_params.μ.u
const DENSITY_THERMO = BENCH_PARAMS.thermo_params
const DENSITY_COS_GRID, DENSITY_COS_W = gauleg(0.0, 1.0, ASR.DEFAULT_ANGLE_NODES)
const ASR_COS_GRID, ASR_COS_W = gauleg(-1.0, 1.0, ASR.DEFAULT_ANGLE_NODES)
const ASR_PHI_GRID, ASR_PHI_W = gauleg(0.0, 2π, ASR.DEFAULT_PHI_NODES)

function register_suite!()
    isdefined(Main, :SUITE) || return nothing

    Main.SUITE["relaxtime"]["number_density_default"] = @benchmarkable ASR.number_density(
        $DENSITY_FLAVOR,
        $DENSITY_MASS,
        $DENSITY_MU,
        $(DENSITY_THERMO.T),
        $(DENSITY_THERMO.Φ),
        $(DENSITY_THERMO.Φbar),
        $(DENSITY_THERMO.ξ);
        p_nodes=$(ASR.DEFAULT_SEMI_INF_NODES),
        angle_nodes=$(ASR.DEFAULT_ANGLE_NODES),
    ) evals=1

    Main.SUITE["relaxtime"]["average_scattering_rate_constant_cache"] = @benchmarkable ASR.average_scattering_rate(
        $BENCH_PROCESS,
        $(BENCH_PARAMS.quark_params),
        $(BENCH_PARAMS.thermo_params),
        $(BENCH_PARAMS.K_coeffs);
        p_nodes=$(ASR.DEFAULT_P_NODES),
        angle_nodes=$(ASR.DEFAULT_ANGLE_NODES),
        phi_nodes=$(ASR.DEFAULT_PHI_NODES),
        cs_cache=$BENCH_CACHE,
        n_sigma_points=4,
    ) evals=1

    Main.SUITE["relaxtime"]["semi_infinite_grid_density_kernel"] = @benchmarkable build_semi_infinite_momentum_grid(
        $(ASR.DEFAULT_SEMI_INF_NODES),
        $(ASR.DEFAULT_SEMI_INF_SCALE),
    ) evals=1

    return nothing
end

function run_report(; density_samples::Int=20, asr_samples::Int=5, kernel_samples::Int=200)
    density_default = @benchmark ASR.number_density(
        $DENSITY_FLAVOR,
        $DENSITY_MASS,
        $DENSITY_MU,
        $(DENSITY_THERMO.T),
        $(DENSITY_THERMO.Φ),
        $(DENSITY_THERMO.Φbar),
        $(DENSITY_THERMO.ξ);
        p_nodes=$(ASR.DEFAULT_SEMI_INF_NODES),
        angle_nodes=$(ASR.DEFAULT_ANGLE_NODES),
    ) samples=density_samples evals=1

    density_reuse_angles = @benchmark ASR.number_density(
        $DENSITY_FLAVOR,
        $DENSITY_MASS,
        $DENSITY_MU,
        $(DENSITY_THERMO.T),
        $(DENSITY_THERMO.Φ),
        $(DENSITY_THERMO.Φbar),
        $(DENSITY_THERMO.ξ);
        p_nodes=$(ASR.DEFAULT_SEMI_INF_NODES),
        angle_nodes=$(ASR.DEFAULT_ANGLE_NODES),
        cos_grid=$DENSITY_COS_GRID,
        cos_w=$DENSITY_COS_W,
    ) samples=density_samples evals=1

    asr_default = @benchmark ASR.average_scattering_rate(
        $BENCH_PROCESS,
        $(BENCH_PARAMS.quark_params),
        $(BENCH_PARAMS.thermo_params),
        $(BENCH_PARAMS.K_coeffs);
        p_nodes=$(ASR.DEFAULT_P_NODES),
        angle_nodes=$(ASR.DEFAULT_ANGLE_NODES),
        phi_nodes=$(ASR.DEFAULT_PHI_NODES),
        cs_cache=$BENCH_CACHE,
        n_sigma_points=4,
    ) samples=asr_samples evals=1

    asr_reuse_angles = @benchmark ASR.average_scattering_rate(
        $BENCH_PROCESS,
        $(BENCH_PARAMS.quark_params),
        $(BENCH_PARAMS.thermo_params),
        $(BENCH_PARAMS.K_coeffs);
        p_nodes=$(ASR.DEFAULT_P_NODES),
        angle_nodes=$(ASR.DEFAULT_ANGLE_NODES),
        phi_nodes=$(ASR.DEFAULT_PHI_NODES),
        cos_grid=$ASR_COS_GRID,
        cos_w=$ASR_COS_W,
        phi_grid=$ASR_PHI_GRID,
        phi_w=$ASR_PHI_W,
        cs_cache=$BENCH_CACHE,
        n_sigma_points=4,
    ) samples=asr_samples evals=1

    density_kernel = @benchmark build_semi_infinite_momentum_grid(
        $(ASR.DEFAULT_SEMI_INF_NODES),
        $(ASR.DEFAULT_SEMI_INF_SCALE),
    ) samples=kernel_samples evals=1

    asr_kernel = @benchmark build_semi_infinite_momentum_grid(
        $(ASR.DEFAULT_P_NODES),
        $(ASR.DEFAULT_SEMI_INF_SCALE),
    ) samples=kernel_samples evals=1

    println("Density + AverageScatteringRate baseline")
    println(repeat("=", 84))
    @printf("process=%s, T=%.1f MeV, muB=%.1f MeV, xi=%.3f\n",
        string(BENCH_PROCESS), 150.0, 800.0, BENCH_PARAMS.thermo_params.ξ)
    @printf("density nodes=(p=%d, angle=%d), asr nodes=(p=%d, angle=%d, phi=%d)\n\n",
        ASR.DEFAULT_SEMI_INF_NODES, ASR.DEFAULT_ANGLE_NODES,
        ASR.DEFAULT_P_NODES, ASR.DEFAULT_ANGLE_NODES, ASR.DEFAULT_PHI_NODES)

    print_summary(summarize_trial("number_density default", density_default))
    print_summary(summarize_trial("number_density reuse angle grid", density_reuse_angles))
    print_summary(summarize_trial("avg_scattering_rate default", asr_default))
    print_summary(summarize_trial("avg_scattering_rate reuse ang grids", asr_reuse_angles))
    print_summary(summarize_trial("semi-inf grid kernel (density)", density_kernel))
    print_summary(summarize_trial("semi-inf grid kernel (asr)", asr_kernel))

    println("\nDecision hints")
    println(repeat("-", 84))
    density_kernel_ratio = median(density_kernel).time / median(density_default).time
    asr_kernel_ratio = median(asr_kernel).time / median(asr_default).time
    @printf("density kernel / number_density median time: %.2f%%\n", 100 * density_kernel_ratio)
    @printf("asr kernel / average_scattering_rate median time: %.2f%%\n", 100 * asr_kernel_ratio)
    return nothing
end

register_suite!()

if abspath(PROGRAM_FILE) == @__FILE__
    density_samples = parse(Int, get(ENV, "BENCH_SAMPLES_DENSITY", "20"))
    asr_samples = parse(Int, get(ENV, "BENCH_SAMPLES_ASR", "5"))
    kernel_samples = parse(Int, get(ENV, "BENCH_SAMPLES_KERNEL", "200"))
    run_report(; density_samples=density_samples, asr_samples=asr_samples, kernel_samples=kernel_samples)
end