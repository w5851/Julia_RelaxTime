# 验证新的默认配置 (HYBRID + n=32)
using Printf
include(joinpath(@__DIR__, "../../src/relaxtime/OneLoopIntegralsAniso.jl"))
include(joinpath(@__DIR__, "../../src/relaxtime/OneLoopIntegrals.jl"))
using QuadGK: quadgk

println("=" ^ 60)
println("Verifying new defaults: STRATEGY_HYBRID + n=32")
println("=" ^ 60)

test_cases = [
    (name="有根(标准)", λ=-1.0, k=0.01, m=0.3, m_prime=0.3),
    (name="有根(大k)", λ=-1.0, k=0.1, m=0.3, m_prime=0.3),
    (name="无根", λ=-0.5, k=0.05, m=0.3, m_prime=0.3),
    (name="无根(正λ)", λ=0.5, k=0.1, m=0.3, m_prime=0.35),
]

ξ = -0.2; T = 0.15; μ = 0.0; Φ = 0.0; Φbar = 0.0

println("\nDefault strategy: ", OneLoopIntegralsCorrection.DEFAULT_STRATEGY)
println("Default n: ", OneLoopIntegralsCorrection.DEFAULT_CLUSTER_N)
println()

all_pass = true
for tc in test_cases
    Emin = tc.m
    Emax = OneLoopIntegrals.energy_cutoff(tc.m)
    integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(:quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
    ref, _ = quadgk(integrand, Emin, Emax; rtol=1e-12)
    
    # 使用默认参数调用
    result = OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
        :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ
    )
    err = abs((result[1] - ref) / ref)
    
    status = err < 1e-4 ? "✓" : "✗"
    global all_pass = all_pass && (err < 1e-4)
    println(@sprintf("%s %-12s: relerr=%.2e %s", status, tc.name, err, err < 1e-4 ? "" : "(FAIL)"))
end

println("\n" * "=" ^ 60)
println(all_pass ? "ALL TESTS PASSED ✓" : "SOME TESTS FAILED ✗")
println("=" ^ 60)
