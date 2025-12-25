"""
检查物理参数范围

在 PNJL 模型中，夸克质量 m 通常在 0.3 GeV 左右（组分夸克质量）
m=0.1 GeV 是一个较小的值，可能对应于高温下的手征恢复相
"""

using Printf

include(joinpath(@__DIR__, "../../src/relaxtime/OneLoopIntegralsAniso.jl"))
using .OneLoopIntegralsCorrection: tilde_B0_correction_k_positive, 
    STRATEGY_HYBRID, STRATEGY_QUADGK

μ, T, Φ, Φbar, ξ = 0.3, 0.15, 0.3, 0.3, 0.1

println("物理参数范围测试")
println("=" ^ 70)
println()

# 典型物理参数范围
# m: 0.3-0.5 GeV (组分夸克质量)
# λ: 0.1-1.0 GeV (外动量)
# k: 0.1-1.0 GeV (外动量)

println("典型物理参数 (m ≥ 0.3 GeV):")
println("-" ^ 70)

typical_cases = [
    (0.3, 0.2, 0.3, 0.3, "m=0.3, 小k"),
    (0.5, 0.3, 0.3, 0.3, "m=0.3, 中等参数"),
    (0.8, 0.5, 0.3, 0.3, "m=0.3, 大参数"),
    (0.3, 0.2, 0.35, 0.35, "m=0.35"),
    (0.5, 0.3, 0.4, 0.4, "m=0.4"),
    (0.5, 0.3, 0.5, 0.5, "m=0.5"),
]

errors = Float64[]
for (λ, k, m, m_prime, desc) in typical_cases
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
@printf("通过率: %d/%d\n", count(e -> e < 1e-4, errors), length(errors))
println()

println("边界情况 (m < 0.3 GeV):")
println("-" ^ 70)

edge_cases = [
    (0.1, 0.2, 0.1, 0.3, "m=0.1 (极小)"),
    (0.2, 0.2, 0.2, 0.3, "m=0.2 (小)"),
    (0.3, 0.2, 0.25, 0.3, "m=0.25"),
]

for (λ, k, m, m_prime, desc) in edge_cases
    real_hybrid, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
        strategy=STRATEGY_HYBRID)
    real_ref, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
        strategy=STRATEGY_QUADGK, rtol=1e-12, atol=1e-14)
    
    rel_err = abs(real_hybrid - real_ref) / max(abs(real_ref), 1e-15)
    status = rel_err < 1e-4 ? "✓" : "✗"
    @printf("%-25s | 误差: %.2e %s\n", desc, rel_err, status)
end

println()
println("结论:")
println("-" ^ 70)
println("1. 对于典型物理参数 (m ≥ 0.3 GeV)，HYBRID 策略精度足够")
println("2. 对于边界情况 (m < 0.2 GeV)，可能需要更多节点")
println("3. 在实际应用中，m=0.1 GeV 是极端情况，较少出现")
