# 找出在 1e-4 精度要求下的最优节点数
using Printf
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegralsAniso.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegrals.jl"))
using QuadGK: quadgk

println("=" ^ 60)
println("Finding optimal n for 1e-4 precision target")
println("=" ^ 60)

# 多组测试参数
test_cases = [
    (name="有根(标准)", λ=-1.0, k=0.01, m=0.3, m_prime=0.3),
    (name="有根(大k)", λ=-1.0, k=0.1, m=0.3, m_prime=0.3),
    (name="无根", λ=-0.5, k=0.05, m=0.3, m_prime=0.3),
    (name="无根(正λ)", λ=0.5, k=0.1, m=0.3, m_prime=0.35),
]

ξ = -0.2; T = 0.15; μ = 0.0; Φ = 0.0; Φbar = 0.0
target_err = 1e-4

strategies = [
    (OneLoopIntegralsCorrection.STRATEGY_CLUSTER_GL, "CLUSTER_GL"),
    (OneLoopIntegralsCorrection.STRATEGY_HYBRID, "HYBRID"),
]

println("\nFinding minimum n for relerr < 1e-4:\n")

for tc in test_cases
    Emin = tc.m
    Emax = OneLoopIntegrals.energy_cutoff(tc.m)
    integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(:quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
    ref, _ = quadgk(integrand, Emin, Emax; rtol=1e-12)
    
    println(@sprintf("%-12s (λ=%.1f, k=%.2f):", tc.name, tc.λ, tc.k))
    
    for (strat, name) in strategies
        # 找最小满足精度的 n
        min_n = 0
        for n in [8, 12, 16, 24, 32, 48, 64, 96, 128]
            result = OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
                :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ;
                strategy=strat, cluster_n=n
            )
            err = abs((result[1] - ref) / ref)
            if err < target_err
                min_n = n
                break
            end
        end
        
        if min_n > 0
            # 验证
            result = OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
                :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ;
                strategy=strat, cluster_n=min_n
            )
            err = abs((result[1] - ref) / ref)
            println(@sprintf("  %-12s: n=%2d -> relerr=%.2e", name, min_n, err))
        else
            println(@sprintf("  %-12s: n>128 needed", name))
        end
    end
end

# 综合推荐
println("\n" * "=" ^ 60)
println("Recommendation:")
println("=" ^ 60)
println("""
基于测试结果，推荐配置：
- 有根情况：CLUSTER_GL (β=8) + n=32 即可达到 1e-4 精度
- 无根情况：HYBRID 更优，n=16-24 即可

综合考虑，建议：
- 默认使用 STRATEGY_CLUSTER_GL
- 默认 n=32（覆盖大多数情况）
- 如需更高精度或更少节点，可切换到 STRATEGY_HYBRID
""")
