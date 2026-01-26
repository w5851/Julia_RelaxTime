#!/usr/bin/env julia
"""
简单验证阈值计算
"""

using Printf

# 使用典型的夸克质量
m_u = 8.134 / 197.33  # MeV -> fm⁻¹
m_s = 206.829 / 197.33  # MeV -> fm⁻¹

println("="^80)
println("阈值验证")
println("="^80)
println()

println("夸克质量:")
@printf("  m_u = %.3f MeV = %.6f fm⁻¹\n", m_u * 197.33, m_u)
@printf("  m_s = %.3f MeV = %.6f fm⁻¹\n", m_s * 197.33, m_s)
println()

# ssbar → uubar
println("过程1: ssbar → uubar")
s_th_init_1 = (m_s + m_s)^2
s_th_final_1 = (m_u + m_u)^2
s_th_1 = max(s_th_init_1, s_th_final_1)

@printf("  初态阈值: (2m_s)² = %.6f fm⁻² = %.1f MeV²\n", s_th_init_1, s_th_init_1 * 197.33^2)
@printf("  末态阈值: (2m_u)² = %.6f fm⁻² = %.1f MeV²\n", s_th_final_1, s_th_final_1 * 197.33^2)
@printf("  总阈值: %.6f fm⁻² = %.1f MeV²\n", s_th_1, s_th_1 * 197.33^2)
println()

# uubar → ssbar
println("过程2: uubar → ssbar")
s_th_init_2 = (m_u + m_u)^2
s_th_final_2 = (m_s + m_s)^2
s_th_2 = max(s_th_init_2, s_th_final_2)

@printf("  初态阈值: (2m_u)² = %.6f fm⁻² = %.1f MeV²\n", s_th_init_2, s_th_init_2 * 197.33^2)
@printf("  末态阈值: (2m_s)² = %.6f fm⁻² = %.1f MeV²\n", s_th_final_2, s_th_final_2 * 197.33^2)
@printf("  总阈值: %.6f fm⁻² = %.1f MeV²\n", s_th_2, s_th_2 * 197.33^2)
println()

println("="^80)
println("结论")
println("="^80)
println()

if abs(s_th_1 - s_th_2) < 1e-10
    println("✓ 阈值相同！")
    @printf("  两个过程的阈值都是 %.1f MeV²\n", s_th_1 * 197.33^2)
    println()
    println("  这是物理正确的：")
    println("  - 对于 A + B ↔ C + D，阈值 = max((m_A+m_B)², (m_C+m_D)²)")
    println("  - ssbar ↔ uubar 的阈值由重粒子决定：(2m_s)²")
else
    println("✗ 阈值不同！（这不应该发生）")
end
println()

println("="^80)
println("之前分析的错误")
println("="^80)
println()

println("我之前错误地说：")
println("  ✗ s_threshold(ssbar→uubar) = 171,113 MeV²")
println("  ✗ s_threshold(uubar→ssbar) = 265 MeV²")
println("  ✗ 比值 = 646")
println()

println("实际上：")
@printf("  ✓ s_threshold(ssbar→uubar) = %.1f MeV²\n", s_th_1 * 197.33^2)
@printf("  ✓ s_threshold(uubar→ssbar) = %.1f MeV²\n", s_th_2 * 197.33^2)
println("  ✓ 比值 = 1.0（完全相同）")
println()

println("这意味着：")
println("  - 两个过程在相同的s值下，都同样接近或远离阈值")
println("  - 不能用\"阈值不同\"来解释3.238倍差异")
println("  - 需要寻找其他原因")
println()

println("="^80)

