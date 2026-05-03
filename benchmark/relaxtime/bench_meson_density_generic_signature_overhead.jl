#!/usr/bin/env julia

# Benchmark: generic signature overhead on meson-density AD key path
#
# Goals:
# - Isolate the cost of widening selected key-path helper signatures from Float64-only
#   to Dual-compatible generic methods.
# - Compare current production generic helpers against local Float64-shadow versions
#   with equivalent formulas.
#
# Run:
#   julia --project=benchmark benchmark/relaxtime/bench_meson_density_generic_signature_overhead.jl

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const BENCH_ENV = normpath(joinpath(@__DIR__, ".."))

pushfirst!(LOAD_PATH, PROJECT_ROOT)
pushfirst!(LOAD_PATH, BENCH_ENV)

using BenchmarkTools
using Printf

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "QuarkDistribution_Aniso.jl"))

using .PNJLQuarkDistributions_Aniso: quark_df_dE, antiquark_df_dE, correction_cos_theta_coefficient

const E0 = 1.2
const P0 = 0.9
const M0 = 0.32
const MU0 = 0.1
const T0 = 0.18
const PHI0 = 0.25
const PHIBAR0 = 0.25
const XI0 = 0.1
const SAMPLES = parse(Int, get(ENV, "BENCH_MD_GENERIC_SAMPLES", "30"))

@inline function _float_exp_clamped(x::Float64)
    return exp(clamp(x, -460.0, 460.0))
end

@inline function quark_df_dE_float(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    β_fm = 1 / T_inv_fm
    exp_term = _float_exp_clamped(-(E_inv_fm - μ_inv_fm) * β_fm)
    exp_term2 = exp_term * exp_term
    exp_term3 = exp_term2 * exp_term

    numerator = Φ * exp_term + 2 * Φbar * exp_term2 + exp_term3
    denominator = 1 + 3 * Φ * exp_term + 3 * Φbar * exp_term2 + exp_term3
    d_numerator = -(Φ * exp_term + 4 * Φbar * exp_term2 + 3 * exp_term3)
    d_denominator = -(3 * Φ * exp_term + 6 * Φbar * exp_term2 + 3 * exp_term3)

    df_dx = d_numerator * denominator - numerator * d_denominator
    df_dx /= denominator ^ 2
    return β_fm * df_dx
end

@inline function antiquark_df_dE_float(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    β_fm = 1 / T_inv_fm
    exp_term = _float_exp_clamped(-(E_inv_fm + μ_inv_fm) * β_fm)
    exp_term2 = exp_term * exp_term
    exp_term3 = exp_term2 * exp_term

    numerator = Φbar * exp_term + 2 * Φ * exp_term2 + exp_term3
    denominator = 1 + 3 * Φbar * exp_term + 3 * Φ * exp_term2 + exp_term3
    d_numerator = -(Φbar * exp_term + 4 * Φ * exp_term2 + 3 * exp_term3)
    d_denominator = -(3 * Φbar * exp_term + 6 * Φ * exp_term2 + 3 * exp_term3)

    df_dx = d_numerator * denominator - numerator * d_denominator
    df_dx /= denominator ^ 2
    return β_fm * df_dx
end

@inline function correction_cos_theta_coefficient_float(
    sign_::Symbol, p_inv_fm::Float64, m_inv_fm::Float64,
    μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64, ξ::Float64)
    E_inv_fm = sqrt(p_inv_fm^2 + m_inv_fm^2)
    coeff = 0.5 * ξ * (p_inv_fm^2) / E_inv_fm
    df_dE = sign_ === :quark ?
        quark_df_dE_float(E_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar) :
        antiquark_df_dE_float(E_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
    return coeff * df_dE
end

function summarize_trial(name::String, trial::BenchmarkTools.Trial)
    estimate = median(trial)
    return (
        name=name,
        time_ns=estimate.time,
        memory_bytes=estimate.memory,
        allocs=estimate.allocs,
    )
end

function print_pair(label::String, generic_summary, float_summary)
    ratio = generic_summary.time_ns / float_summary.time_ns
    @printf("%-28s generic=%10.1f ns  float=%10.1f ns  ratio=%6.3f  alloc=%d/%d\n",
        label,
        generic_summary.time_ns,
        float_summary.time_ns,
        ratio,
        generic_summary.allocs,
        float_summary.allocs,
    )
end

function run_report()
    # warmup
    quark_df_dE(E0, MU0, T0, PHI0, PHIBAR0)
    quark_df_dE_float(E0, MU0, T0, PHI0, PHIBAR0)
    antiquark_df_dE(E0, MU0, T0, PHI0, PHIBAR0)
    antiquark_df_dE_float(E0, MU0, T0, PHI0, PHIBAR0)
    correction_cos_theta_coefficient(:quark, P0, M0, MU0, T0, PHI0, PHIBAR0, XI0)
    correction_cos_theta_coefficient_float(:quark, P0, M0, MU0, T0, PHI0, PHIBAR0, XI0)

    q_generic = @benchmark quark_df_dE($E0, $MU0, $T0, $PHI0, $PHIBAR0) samples=SAMPLES evals=1
    q_float = @benchmark quark_df_dE_float($E0, $MU0, $T0, $PHI0, $PHIBAR0) samples=SAMPLES evals=1

    aq_generic = @benchmark antiquark_df_dE($E0, $MU0, $T0, $PHI0, $PHIBAR0) samples=SAMPLES evals=1
    aq_float = @benchmark antiquark_df_dE_float($E0, $MU0, $T0, $PHI0, $PHIBAR0) samples=SAMPLES evals=1

    corr_generic = @benchmark correction_cos_theta_coefficient(:quark, $P0, $M0, $MU0, $T0, $PHI0, $PHIBAR0, $XI0) samples=SAMPLES evals=1
    corr_float = @benchmark correction_cos_theta_coefficient_float(:quark, $P0, $M0, $MU0, $T0, $PHI0, $PHIBAR0, $XI0) samples=SAMPLES evals=1

    println("Meson-density generic signature overhead benchmark")
    println(repeat("=", 88))
    println("Compare current generic helpers vs local Float64-shadow formulas")
    println()

    print_pair("quark_df_dE", summarize_trial("quark_df_dE_generic", q_generic), summarize_trial("quark_df_dE_float", q_float))
    print_pair("antiquark_df_dE", summarize_trial("antiquark_df_dE_generic", aq_generic), summarize_trial("antiquark_df_dE_float", aq_float))
    print_pair("correction_coeff", summarize_trial("correction_generic", corr_generic), summarize_trial("correction_float", corr_float))
    return nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_report()
end
