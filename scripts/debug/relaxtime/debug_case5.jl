"""
调试 Case 5
"""

using Printf

include("../../../src/relaxtime/OneLoopIntegralsAniso.jl")
using .OneLoopIntegralsCorrection: tilde_B0_correction_k_positive, 
    STRATEGY_HYBRID, STRATEGY_QUADGK, STRATEGY_CLUSTER_GL

# Case 5: λ=0.1, k=0.2, m=0.1, m'=0.3
λ, k, m, m_prime = 0.1, 0.2, 0.1, 0.3
μ, T, Φ, Φbar, ξ = 0.3, 0.15, 0.3, 0.3, 0.1

println("Case 5: λ=$λ, k=$k, m=$m, m'=$m_prime")
println("=" ^ 60)

# 诊断信息
_, _, diag = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=STRATEGY_HYBRID, diagnostics=true)

println("根数: $(diag.n_roots)")
println("区间: $(diag.intervals)")

for (i, (a, b)) in enumerate(diag.intervals)
    len = b - a
    @printf("  区间 %d: [%.6f, %.6f], 长度=%.6f\n", i, a, b, len)
end

# 参考值
real_ref, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=STRATEGY_QUADGK, rtol=1e-12, atol=1e-14)

println("\n参考值: $real_ref")

# 不同策略比较
strategies = [
    (STRATEGY_HYBRID, "HYBRID"),
    (STRATEGY_QUADGK, "QuadGK (默认)"),
    (STRATEGY_CLUSTER_GL, "CLUSTER_GL (n=32)"),
]

println("\n策略比较:")
println("-" ^ 60)

for (strat, name) in strategies
    if strat == STRATEGY_CLUSTER_GL
        real_val, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=strat, cluster_n=32)
    else
        real_val, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=strat)
    end
    rel_err = abs(real_val - real_ref) / abs(real_ref)
    @printf("%-20s: %.15e (误差: %.2e)\n", name, real_val, rel_err)
end

# 测试更多节点
println("\nCLUSTER_GL 不同节点数:")
println("-" ^ 60)

for n in [32, 48, 64, 96, 128]
    real_val, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
        strategy=STRATEGY_CLUSTER_GL, cluster_n=n)
    rel_err = abs(real_val - real_ref) / abs(real_ref)
    @printf("n=%3d: %.15e (误差: %.2e)\n", n, real_val, rel_err)
end

println("\n注意: 这个用例的 m=0.1 很小，导致积分区间左端点 E=m=0.1 处")
println("p = sqrt(E²-m²) 的奇异行为更加剧烈")
