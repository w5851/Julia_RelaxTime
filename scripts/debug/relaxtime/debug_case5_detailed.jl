"""
详细调试 Case 5 - 检查 HYBRID 的节点分配
"""

using Printf

include("../../../src/relaxtime/OneLoopIntegralsAniso.jl")
using .OneLoopIntegralsCorrection: tilde_B0_correction_k_positive, 
    STRATEGY_HYBRID, STRATEGY_QUADGK

# Case 5: λ=0.1, k=0.2, m=0.1, m'=0.3
λ, k, m, m_prime = 0.1, 0.2, 0.1, 0.3
μ, T, Φ, Φbar, ξ = 0.3, 0.15, 0.3, 0.3, 0.1

println("Case 5 详细分析")
println("=" ^ 60)

# 诊断信息
_, _, diag = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=STRATEGY_HYBRID, diagnostics=true)

println("根数: $(diag.n_roots)")
println("区间数: $(diag.n_intervals)")
println()

# 当前节点分配（根据代码逻辑）
# - 第一个区间 (SING_RIGHT): 32 节点
# - 最后一个区间 (SING_LEFT): 16 节点
# 总计: 32 + 16 = 48 节点

println("当前节点分配:")
println("  区间 1 (SING_RIGHT): 32 节点")
println("  区间 2 (SING_LEFT):  16 节点")
println("  总计: 48 节点")
println()

# 参考值
real_ref, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=STRATEGY_QUADGK, rtol=1e-12, atol=1e-14)

real_hybrid, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=STRATEGY_HYBRID)

println("参考值:  $real_ref")
println("HYBRID:  $real_hybrid")
println("相对误差: $(abs(real_hybrid - real_ref) / abs(real_ref))")
println()

# 问题分析
println("问题分析:")
println("-" ^ 60)
println("1. m=0.1 很小，导致 E=m 处的 p→0 奇异行为更剧烈")
println("2. 第一个区间 [0.1, 0.201] 长度约 0.1，但包含了")
println("   从 E=m 到根位置的整个范围")
println("3. 被积函数在 E→m 时有 (E-m)^(-1/2) 的奇异行为")
println()

# 检查是否是第一个区间的问题
# 通过比较不同策略在第一个区间的表现
println("建议:")
println("1. 对于 m 很小的情况，可能需要更多节点")
println("2. 或者调整 power_left 变换的 alpha 参数")
println("3. 当前 alpha=0.35，可以尝试更小的值（更强的聚簇）")
