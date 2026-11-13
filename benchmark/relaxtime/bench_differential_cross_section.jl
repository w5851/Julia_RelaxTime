#!/usr/bin/env julia
# Run: julia --project=. benchmark/relaxtime/bench_differential_cross_section.jl

using BenchmarkTools

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxTime.jl"))

function run_report(; samples::Int=100)
    samples > 0 || throw(ArgumentError("samples must be positive"))
    T = 150.0 / 197.327
    m_u, m_s = 300.0 / 197.327, 500.0 / 197.327
    nodes, weights = GaussLegendre.gauleg(0.0, 20.0, 64)
    A_u = OneLoopIntegrals.A(m_u, 0.0, T, 0.5, 0.5, nodes, weights)
    A_s = OneLoopIntegrals.A(m_s, 0.0, T, 0.5, 0.5, nodes, weights)
    K = EffectiveCouplings.calculate_effective_couplings(
        Constants_PNJL.G_fm2, Constants_PNJL.K_fm5,
        EffectiveCouplings.calculate_G_from_A(A_u, m_u),
        EffectiveCouplings.calculate_G_from_A(A_s, m_s),
    )
    quark = (m=(u=m_u, d=m_u, s=m_s), μ=(u=0.0, d=0.0, s=0.0), A=(u=A_u, d=A_u, s=A_s))
    thermo = (T=T, Φ=0.5, Φbar=0.5, ξ=0.0)
    s, t = 31.0, -2.0
    process = :uu_to_uu
    variables = ScatteringAmplitude.calculate_mandelstam_variables(s, t, 4m_u^2 - s - t, m_u, m_u, m_u, m_u)
    amplitude = ScatteringAmplitude.scattering_amplitude_squared(process, s, t, quark, thermo, K)
    suite = BenchmarkGroup()
    suite["differential_cross_section"] = @benchmarkable DifferentialCrossSection.differential_cross_section(
        $(Ref(variables.s_12_plus))[], $(Ref(variables.s_12_minus))[], $(Ref(amplitude))[]) samples=samples evals=1000
    suite["kinematic_threshold"] = @benchmarkable DifferentialCrossSection.check_kinematic_threshold(
        $(Ref(s))[], $(Ref(m_u))[], $(Ref(m_u))[]; warn_close=false) samples=samples evals=1000
    suite["t_bounds"] = @benchmarkable TotalCrossSection.calculate_t_bounds(
        $(Ref(s))[], $(Ref(m_u))[], $(Ref(m_u))[], $(Ref(m_u))[], $(Ref(m_u))[]) samples=samples evals=1000
    suite["amplitude_and_cross_section"] = @benchmarkable DifferentialCrossSection.differential_cross_section(
        $(variables.s_12_plus), $(variables.s_12_minus),
        ScatteringAmplitude.scattering_amplitude_squared($process, $s, $t, $quark, $thermo, $K)) samples=samples evals=1
    result = run(suite)
    show(stdout, MIME("text/plain"), result)
    println()
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_report(; samples=parse(Int, get(ENV, "BENCH_SAMPLES", "100")))
end
