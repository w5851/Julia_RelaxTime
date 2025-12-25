"""
调试精度问题
"""

using Printf

include("../../../src/relaxtime/OneLoopIntegralsAniso.jl")
using .OneLoopIntegralsCorrection: tilde_B0_correction_k_positive, 
    STRATEGY_HYBRID, STRATEGY_QUADGK, STRATEGY_CLUSTER_GL

# 问题用例
λ, k, m, m_prime = 0.1, 0.3, 0.3, 0.5
μ, T, Φ, Φbar, ξ = 0.3, 0.15, 0.3, 0.3, 0.1

println("问题用例: λ=$λ, k=$k, m=$m, m'=$m_prime")
println("=" ^ 60)

# 获取诊断信息
_, _, diag = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=STRATEGY_HYBRID, diagnostics=true)

println("根数: $(diag.n_roots)")
println("区间数: $(diag.n_intervals)")
println("根位置: $(diag.roots)")
println("区间: $(diag.intervals)")
println()

# 比较不同策略
real_hybrid, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=STRATEGY_HYBRID)

real_quadgk_default, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=STRATEGY_QUADGK)

real_quadgk_high, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=STRATEGY_QUADGK, rtol=1e-10, atol=1e-12)

real_cluster_32, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=STRATEGY_CLUSTER_GL, cluster_n=32)

real_cluster_64, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=STRATEGY_CLUSTER_GL, cluster_n=64)

println("结果比较:")
println("-" ^ 60)
@printf("HYBRID (自适应):     %.15e\n", real_hybrid)
@printf("QuadGK (默认):       %.15e\n", real_quadgk_default)
@printf("QuadGK (高精度):     %.15e\n", real_quadgk_high)
@printf("CLUSTER_GL (n=32):   %.15e\n", real_cluster_32)
@printf("CLUSTER_GL (n=64):   %.15e\n", real_cluster_64)
println()

# 以高精度 QuadGK 为参考
ref = real_quadgk_high
println("相对误差 (参考: QuadGK 高精度):")
println("-" ^ 60)
@printf("HYBRID:              %.2e\n", abs(real_hybrid - ref) / abs(ref))
@printf("QuadGK (默认):       %.2e\n", abs(real_quadgk_default - ref) / abs(ref))
@printf("CLUSTER_GL (n=32):   %.2e\n", abs(real_cluster_32 - ref) / abs(ref))
@printf("CLUSTER_GL (n=64):   %.2e\n", abs(real_cluster_64 - ref) / abs(ref))
println()

# 检查 QuadGK 默认精度是否足够
println("注意: QuadGK 默认精度的误差也较大，说明这是一个困难的积分")
println("可能需要增加有根情况下的节点数")
