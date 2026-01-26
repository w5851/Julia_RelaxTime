#!/usr/bin/env julia
"""
测试：使用 Fortran 的 m, Φ 计算 A 函数

目的：验证 Julia 的 A 函数实现是否与 Fortran 一致

测试条件：T=300 MeV, μ=2 MeV
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
println("使用 Fortran 的 m, Φ 计算 A 函数")
println("="^80)
println()

# 物理参数
T_MeV = 300.0
μ_MeV = 2.0

T = T_MeV / ħc_MeV_fm
μ = μ_MeV / ħc_MeV_fm

@printf("物理参数:\n")
@printf("  T = %.1f MeV = %.6f fm⁻¹\n", T_MeV, T)
@printf("  μ = %.1f MeV = %.6f fm⁻¹\n", μ_MeV, μ)
println()

# Fortran 的解
m_u_fortran = 0.040510  # fm⁻¹
m_d_fortran = 0.040510  # fm⁻¹
m_s_fortran = 1.032296  # fm⁻¹
Φ_fortran = 0.99999994
Φbar_fortran = 0.99999994

println("Fortran 的平衡态:")
@printf("  m_u     = %.6f fm⁻¹ = %.2f MeV\n", m_u_fortran, m_u_fortran * ħc_MeV_fm)
@printf("  m_d     = %.6f fm⁻¹ = %.2f MeV\n", m_d_fortran, m_d_fortran * ħc_MeV_fm)
@printf("  m_s     = %.6f fm⁻¹ = %.2f MeV\n", m_s_fortran, m_s_fortran * ħc_MeV_fm)
@printf("  Φ       = %.8f\n", Φ_fortran)
@printf("  Φ̄       = %.8f\n", Φbar_fortran)
println()

# 生成积分节点（与 Fortran 一致：128 个节点）
nodes_p, weights_p = gauleg(0.0, 20.0, 128)

println("="^80)
println("计算 A 函数")
println("="^80)
println()

# 计算 A 函数
A_u = OneLoopIntegrals.A(m_u_fortran, μ, T, Φ_fortran, Φbar_fortran, nodes_p, weights_p)
A_d = OneLoopIntegrals.A(m_d_fortran, μ, T, Φ_fortran, Φbar_fortran, nodes_p, weights_p)
A_s = OneLoopIntegrals.A(m_s_fortran, μ, T, Φ_fortran, Φbar_fortran, nodes_p, weights_p)

println("Julia 计算的 A 函数:")
@printf("  A_u     = %.6f fm⁻²\n", A_u)
@printf("  A_d     = %.6f fm⁻²\n", A_d)
@printf("  A_s     = %.6f fm⁻²\n", A_s)
println()

# Fortran 的 A 值
A_u_fortran = -4.972144
A_d_fortran = -4.972144
A_s_fortran = -5.165541

println("Fortran 的 A 函数:")
@printf("  A_u     = %.6f fm⁻²\n", A_u_fortran)
@printf("  A_d     = %.6f fm⁻²\n", A_d_fortran)
@printf("  A_s     = %.6f fm⁻²\n", A_s_fortran)
println()

# 对比
println("="^80)
println("对比分析")
println("="^80)
println()

diff_u = abs(A_u - A_u_fortran) / abs(A_u_fortran) * 100
diff_d = abs(A_d - A_d_fortran) / abs(A_d_fortran) * 100
diff_s = abs(A_s - A_s_fortran) / abs(A_s_fortran) * 100

@printf("相对差异:\n")
@printf("  ΔA_u    = %.2f%%\n", diff_u)
@printf("  ΔA_d    = %.2f%%\n", diff_d)
@printf("  ΔA_s    = %.2f%%\n", diff_s)
println()

# 计算能隙方程残差
println("="^80)
println("检查能隙方程残差")
println("="^80)
println()

residual_u = m_u_fortran - m_ud0_inv_fm + 2.0 * G_fm2 * m_u_fortran * A_u - 
             4.0 * K_fm5 * m_u_fortran * A_u * A_d * A_s
residual_d = m_d_fortran - m_ud0_inv_fm + 2.0 * G_fm2 * m_d_fortran * A_d - 
             4.0 * K_fm5 * m_d_fortran * A_u * A_d * A_s
residual_s = m_s_fortran - m_s0_inv_fm + 2.0 * G_fm2 * m_s_fortran * A_s - 
             4.0 * K_fm5 * m_s_fortran * A_u * A_d * A_s

println("能隙方程残差（使用 Julia 的 A）:")
@printf("  F_u     = %.2e\n", residual_u)
@printf("  F_d     = %.2e\n", residual_d)
@printf("  F_s     = %.2e\n", residual_s)
println()

# 使用 Fortran 的 A 计算残差
residual_u_f = m_u_fortran - m_ud0_inv_fm + 2.0 * G_fm2 * m_u_fortran * A_u_fortran - 
               4.0 * K_fm5 * m_u_fortran * A_u_fortran * A_d_fortran * A_s_fortran
residual_d_f = m_d_fortran - m_ud0_inv_fm + 2.0 * G_fm2 * m_d_fortran * A_d_fortran - 
               4.0 * K_fm5 * m_d_fortran * A_u_fortran * A_d_fortran * A_s_fortran
residual_s_f = m_s_fortran - m_s0_inv_fm + 2.0 * G_fm2 * m_s_fortran * A_s_fortran - 
               4.0 * K_fm5 * m_s_fortran * A_u_fortran * A_d_fortran * A_s_fortran

println("能隙方程残差（使用 Fortran 的 A）:")
@printf("  F_u     = %.2e\n", residual_u_f)
@printf("  F_d     = %.2e\n", residual_d_f)
@printf("  F_s     = %.2e\n", residual_s_f)
println()

# 总结
println("="^80)
println("总结")
println("="^80)
println()

if diff_u < 1.0 && diff_s < 1.0
    println("✅ Julia 的 A 函数与 Fortran 一致（差异 < 1%）")
    println("   → A 函数的实现是正确的")
    println("   → Fortran 的解确实满足能隙方程")
    println()
    println("   结论：问题不在 A 函数，而在平衡态求解器")
    println("   → Julia 的求解器收敛到了不同的局部解")
    println("   → 需要检查初始值策略或求解器算法")
else
    println("❌ Julia 的 A 函数与 Fortran 不一致（差异 > 1%）")
    println("   → A 函数的实现可能不同")
    println()
    println("   可能的原因：")
    println("   1. 积分方法不同（节点数、范围）")
    println("   2. 分布函数的实现不同")
    println("   3. 常数项的计算不同")
    println("   4. 数值精度问题")
end
println()

if abs(residual_u_f) < 1e-6 && abs(residual_s_f) < 1e-6
    println("✅ Fortran 的解满足能隙方程（残差 < 1e-6）")
    println("   → Fortran 的解是正确的")
else
    println("⚠️  Fortran 的解不完全满足能隙方程")
    println("   → 可能是收敛标准不够严格")
end
println()

println("="^80)
