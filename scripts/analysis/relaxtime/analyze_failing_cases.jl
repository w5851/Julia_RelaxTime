"""
分析精度不足的用例
"""

using Printf

include("../../../src/relaxtime/OneLoopIntegralsAniso.jl")
using .OneLoopIntegralsCorrection: tilde_B0_correction_k_positive, 
    STRATEGY_HYBRID, STRATEGY_QUADGK

μ, T, Φ, Φbar, ξ = 0.3, 0.15, 0.3, 0.3, 0.1

# 问题用例
failing_cases = [
    (0.1, 0.3, 0.3, 0.5, "Case 2"),
    (0.1, 0.2, 0.1, 0.3, "Case 5"),
]

println("分析精度不足的用例")
println("=" ^ 70)

for (λ, k, m, m_prime, desc) in failing_cases
    println("\n$desc: λ=$λ, k=$k, m=$m, m'=$m_prime")
    println("-" ^ 70)
    
    _, _, diag = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
        strategy=STRATEGY_HYBRID, diagnostics=true)
    
    println("根数: $(diag.n_roots)")
    println("区间: $(diag.intervals)")
    
    # 计算区间长度
    for (i, (a, b)) in enumerate(diag.intervals)
        len = b - a
        @printf("  区间 %d: [%.6f, %.6f], 长度=%.6f\n", i, a, b, len)
    end
    
    # 检查根位置相对于区间的位置
    if !isempty(diag.roots)
        root = diag.roots[1]
        Emin = m
        Emax = 3.0670014837478465  # energy_cutoff(m)
        rel_pos = (root - Emin) / (Emax - Emin)
        @printf("根位置: %.6f (相对位置: %.2f%%)\n", root, rel_pos * 100)
    end
    
    # 参考值
    real_ref, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
        strategy=STRATEGY_QUADGK, rtol=1e-12, atol=1e-14)
    
    real_hybrid, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
        strategy=STRATEGY_HYBRID)
    
    rel_err = abs(real_hybrid - real_ref) / abs(real_ref)
    @printf("相对误差: %.2e\n", rel_err)
end

println("\n" * "=" ^ 70)
println("分析结论:")
println("-" ^ 70)
println("问题用例的共同特点:")
println("1. 根位置靠近积分区间左端点 (E=m)")
println("2. 第一个区间很短，但被积函数在此区间变化剧烈")
println("3. 需要在第一个区间使用更多节点")
println()
println("建议: 对于有根情况，将第一个区间的节点数从 16 增加到 32")
