#!/usr/bin/env julia
"""
验证积分范围对 A 函数的影响

测试不同积分范围下的 A 函数值，验证这是否是 31% 差异的原因
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))

using .Constants_PNJL
using .GaussLegendre
using .OneLoopIntegrals
using Printf

println("="^80)
println("验证积分范围对 A 函数的影响")
println("="^80)
println()

# 物理参数
T_MeV = 300.0
μ_MeV = 2.0
T = T_MeV / ħc_MeV_fm
μ = μ_MeV / ħc_MeV_fm

# Fortran 的平衡态
m_u_fortran = 0.040510  # fm⁻¹
m_s_fortran = 1.032296  # fm⁻¹
Φ_fortran = 0.99999994
Φbar_fortran = 0.99999994

@printf("物理参数:\n")
@printf("  T = %.1f MeV\n", T_MeV)
@printf("  μ = %.1f MeV\n", μ_MeV)
@printf("  m_u = %.6f fm⁻¹\n", m_u_fortran)
@printf("  m_s = %.6f fm⁻¹\n", m_s_fortran)
println()

# Fortran 的 A 值
A_u_fortran = -4.972144
A_s_fortran = -5.165541

println("="^80)
println("测试不同积分范围")
println("="^80)
println()

# 测试不同的积分范围
ranges = [10.0, 15.0, 20.0, 25.0, 30.0]

println("A_u 随积分范围的变化:")
println("-"^80)
@printf("%-15s %-20s %-20s %-15s\n", "范围 (fm⁻¹)", "A_u (Julia)", "A_u (Fortran)", "差异 (%)")
println("-"^80)

for p_max in ranges
    nodes_p, weights_p = gauleg(0.0, p_max, 128)
    A_u = OneLoopIntegrals.A(m_u_fortran, μ, T, Φ_fortran, Φbar_fortran, nodes_p, weights_p)
    diff = abs(A_u - A_u_fortran) / abs(A_u_fortran) * 100
    @printf("%-15.1f %-20.6f %-20.6f %-15.2f\n", p_max, A_u, A_u_fortran, diff)
end
println()

println("A_s 随积分范围的变化:")
println("-"^80)
@printf("%-15s %-20s %-20s %-15s\n", "范围 (fm⁻¹)", "A_s (Julia)", "A_s (Fortran)", "差异 (%)")
println("-"^80)

for p_max in ranges
    nodes_p, weights_p = gauleg(0.0, p_max, 128)
    A_s = OneLoopIntegrals.A(m_s_fortran, μ, T, Φ_fortran, Φbar_fortran, nodes_p, weights_p)
    diff = abs(A_s - A_s_fortran) / abs(A_s_fortran) * 100
    @printf("%-15.1f %-20.6f %-20.6f %-15.2f\n", p_max, A_s, A_s_fortran, diff)
end
println()

# 重点测试 Fortran 的范围
println("="^80)
println("使用 Fortran 的积分范围 [0, 15 fm⁻¹]")
println("="^80)
println()

nodes_p_15, weights_p_15 = gauleg(0.0, 15.0, 128)
A_u_15 = OneLoopIntegrals.A(m_u_fortran, μ, T, Φ_fortran, Φbar_fortran, nodes_p_15, weights_p_15)
A_s_15 = OneLoopIntegrals.A(m_s_fortran, μ, T, Φ_fortran, Φbar_fortran, nodes_p_15, weights_p_15)

println("Julia 结果 (15 fm⁻¹):")
@printf("  A_u = %.6f fm⁻²\n", A_u_15)
@printf("  A_s = %.6f fm⁻²\n", A_s_15)
println()

println("Fortran 结果:")
@printf("  A_u = %.6f fm⁻²\n", A_u_fortran)
@printf("  A_s = %.6f fm⁻²\n", A_s_fortran)
println()

println("差异:")
diff_u_15 = abs(A_u_15 - A_u_fortran) / abs(A_u_fortran) * 100
diff_s_15 = abs(A_s_15 - A_s_fortran) / abs(A_s_fortran) * 100
@printf("  ΔA_u = %.2f%%\n", diff_u_15)
@printf("  ΔA_s = %.2f%%\n", diff_s_15)
println()

# 总结
println("="^80)
println("总结")
println("="^80)
println()

if diff_u_15 < 5.0 && diff_s_15 < 5.0
    println("✅ 使用 Fortran 的积分范围 [0, 15 fm⁻¹]，差异显著减小（< 5%）")
    println("   → 积分范围是主要差异来源")
    println()
    println("   但仍有 ~", round(diff_u_15, digits=1), "% 的差异，可能来自：")
    println("   1. 常数项的计算方法（Fortran 数值 vs Julia 解析）")
    println("   2. 分布函数的实现细节")
    println("   3. 数值精度")
elseif diff_u_15 < 1.0 && diff_s_15 < 1.0
    println("✅✅ 使用 Fortran 的积分范围 [0, 15 fm⁻¹]，差异几乎消失（< 1%）")
    println("   → 积分范围是唯一的主要差异来源")
    println("   → Julia 的实现是正确的")
    println()
    println("   建议：")
    println("   - 统一 Fortran 和 Julia 的积分范围")
    println("   - 或都使用半无穷积分")
    println("   - 或根据温度自适应选择范围")
else
    println("❌ 即使使用 Fortran 的积分范围，差异仍然很大（> 5%）")
    println("   → 积分范围不是唯一原因")
    println()
    println("   需要进一步检查：")
    println("   1. 常数项的计算")
    println("   2. 分布函数的实现")
    println("   3. 角度积分的处理")
end
println()

# 检查收敛性
println("="^80)
println("检查积分收敛性")
println("="^80)
println()

println("A_u 的收敛性:")
for i in 1:length(ranges)-1
    nodes_1, weights_1 = gauleg(0.0, ranges[i], 128)
    nodes_2, weights_2 = gauleg(0.0, ranges[i+1], 128)
    A_1 = OneLoopIntegrals.A(m_u_fortran, μ, T, Φ_fortran, Φbar_fortran, nodes_1, weights_1)
    A_2 = OneLoopIntegrals.A(m_u_fortran, μ, T, Φ_fortran, Φbar_fortran, nodes_2, weights_2)
    change = abs(A_2 - A_1) / abs(A_1) * 100
    @printf("  [%.1f → %.1f] fm⁻¹: 变化 %.2f%%\n", ranges[i], ranges[i+1], change)
end
println()

if abs(A_u_15 - A_u_fortran) / abs(A_u_fortran) < 0.01
    println("✅ 在 15 fm⁻¹ 处已经收敛")
else
    println("⚠️  在 15 fm⁻¹ 处尚未完全收敛")
    println("   建议使用更大的积分范围或半无穷积分")
end
println()

println("="^80)
