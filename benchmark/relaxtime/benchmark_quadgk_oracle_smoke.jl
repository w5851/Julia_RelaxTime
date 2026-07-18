# Isolated QuadGK oracle smoke for the current fixed-hybrid B0 implementation.
#
# Required environment setup is documented in benchmark/README.md.

using QuadGK: quadgk

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxTime.jl"))

const OLI = Main.OneLoopIntegrals
const OLIC = Main.OneLoopIntegralsCorrection

const CASES = [
    (name="root-sensitive", λ=-1.0, k=0.1, m=0.3, m_prime=0.3),
    (name="smooth", λ=-0.5, k=0.05, m=0.3, m_prime=0.3),
]

function main()
    ξ = -0.2
    T = 0.15
    μ = 0.0
    Φ = 0.0
    Φbar = 0.0

    for case in CASES
        Emin = case.m
        Emax = OLI.energy_cutoff(case.m)
        integrand(E) = OLIC.real_integrand_k_positive(
            :quark, case.λ, case.k, case.m, case.m_prime, E, ξ, T, μ, Φ, Φbar,
        )
        roots = OLIC.find_roots_AB(case.λ, case.k, case.m, case.m_prime, Emin, Emax)
        points = Float64[Emin; roots; Emax]
        oracle, error_estimate = quadgk(integrand, points...; rtol=1e-9, atol=1e-11)
        hybrid, _ = OLIC.tilde_B0_correction_k_positive(
            :quark, case.λ, case.k, case.m, case.m_prime, μ, T, Φ, Φbar, ξ,
        )
        relative_difference = abs(hybrid - oracle) / max(abs(oracle), eps(Float64))
        println((;
            case=case.name,
            roots=length(roots),
            hybrid,
            oracle,
            error_estimate,
            relative_difference,
        ))
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
