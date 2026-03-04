"""
诊断脚本：测试 OneLoopIntegralsAniso 的不同节点分配/求积策略在若干有根用例上的误差表现。

注意：这是诊断/分析脚本，不属于单元/集成测试套件。
运行方式示例：
  julia --project=. tests/analysis/relaxtime_diagnostics/node_allocation_diagnostic.jl
"""

using Printf

# 从仓库根目录定位 src
const _ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(_ROOT, "src", "relaxtime", "OneLoopIntegralsAniso.jl"))

using .OneLoopIntegralsCorrection: tilde_B0_correction_k_positive,
    STRATEGY_HYBRID, STRATEGY_QUADGK

# 问题用例
λ, k, m, m_prime = 0.1, 0.3, 0.3, 0.5
μ, T, Φ, Φbar, ξ = 0.3, 0.15, 0.3, 0.3, 0.1

# 参考值
real_ref, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=STRATEGY_QUADGK, rtol=1e-12, atol=1e-14)

println("参考值 (QuadGK 高精度): $real_ref")
println()

# 当前 HYBRID 结果
real_hybrid, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=STRATEGY_HYBRID)

println("当前 HYBRID (默认节点): $real_hybrid")
println("相对误差: $(abs(real_hybrid - real_ref) / max(abs(real_ref), 1e-15))")
println()

# 测试更多用例
println("=" ^ 70)
println("测试更多有根用例:")
println("-" ^ 70)

test_cases = [
    (0.1, 0.2, 0.3, 0.3, "Case 1"),
    (0.1, 0.3, 0.3, 0.5, "Case 2 (问题用例)"),
    (0.2, 0.3, 0.3, 0.3, "Case 3"),
    (0.3, 0.5, 0.3, 0.3, "Case 4"),
    (0.1, 0.2, 0.1, 0.3, "Case 5"),
    (0.2, 0.4, 0.3, 0.5, "Case 6"),
]

errors = Float64[]
for (λ, k, m, m_prime, desc) in test_cases
    real_hybrid, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
        strategy=STRATEGY_HYBRID)
    real_ref, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
        strategy=STRATEGY_QUADGK, rtol=1e-12, atol=1e-14)

    rel_err = abs(real_hybrid - real_ref) / max(abs(real_ref), 1e-15)
    push!(errors, rel_err)

    status = rel_err < 1e-4 ? "✓" : "✗"
    @printf("%-25s | 误差: %.2e %s\n", desc, rel_err, status)
end

println("-" ^ 70)
@printf("最大误差: %.2e\n", maximum(errors))
@printf("平均误差: %.2e\n", sum(errors) / length(errors))
println()

pass_count = count(e -> e < 1e-4, errors)
@printf("通过率: %d/%d (%.1f%%)\n", pass_count, length(errors), 100 * pass_count / length(errors))
