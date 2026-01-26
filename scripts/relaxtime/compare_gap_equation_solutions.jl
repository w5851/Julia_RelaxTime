#!/usr/bin/env julia
"""
对比 Fortran 和 Julia 的能隙方程解（序参量）

从 Fortran 调试输出:
  m_u = 0.04051 fm⁻¹
  m_s = 1.0323 fm⁻¹
  
需要检查:
1. 这些质量是如何从能隙方程求解得到的
2. Julia 使用相同的初始条件能否得到相同的解
3. 序参量 φ_u, φ_d, φ_s 是否一致
4. Polyakov 环 Φ, Φbar 是否一致
"""

using Printf

# 包含必要的模块
include("../../src/Constants_PNJL.jl")
using .Constants_PNJL

# 测试参数
const ħc = 0.1973269804  # GeV·fm
T_MeV = 300.0
μ_MeV = 2.0
T = T_MeV / (ħc * 1000.0)
μ = μ_MeV / (ħc * 1000.0)

println("=" ^ 80)
println("能隙方程解对比")
println("=" ^ 80)
println()
println("测试条件:")
@printf("  T = %.6f fm⁻¹ (%.1f MeV)\n", T, T_MeV)
@printf("  μ = %.6f fm⁻¹ (%.1f MeV)\n", μ, μ_MeV)
println()

println("=" ^ 80)
println("Fortran 的解 (从调试输出)")
println("=" ^ 80)
println()

# 从 Fortran 调试输出
m_u_fortran = 0.04051  # fm⁻¹
m_s_fortran = 1.0323   # fm⁻¹
Phi_fortran = 0.99999994
Phibar_fortran = 0.99999994

@printf("质量:\n")
@printf("  m_u = %.5f fm⁻¹ (%.2f MeV)\n", m_u_fortran, m_u_fortran * ħc * 1000.0)
@printf("  m_s = %.5f fm⁻¹ (%.2f MeV)\n", m_s_fortran, m_s_fortran * ħc * 1000.0)
println()

@printf("Polyakov 环:\n")
@printf("  Φ = %.8f\n", Phi_fortran)
@printf("  Φbar = %.8f\n", Phibar_fortran)
println()

# 从质量反推序参量
# m = m_0 - 4G⟨ψ̄ψ⟩ + 2K⟨ψ̄ψ⟩_other
# 对于 u 夸克: m_u = m_0 - 4G φ_u + 2K φ_d φ_s
# 对于 s 夸克: m_s = m_0 - 4G φ_s + 2K φ_u φ_d

# 从 Constants_PNJL 获取参数
m_0 = m0_u_inv_fm  # 裸质量
G = G_inv_fm2      # 四费米耦合常数
K = K_inv_fm5      # 六费米耦合常数

@printf("物理参数:\n")
@printf("  m_0 = %.5f fm⁻¹ (%.2f MeV)\n", m_0, m_0 * ħc * 1000.0)
@printf("  G = %.6f fm² (%.2f MeV⁻²)\n", G, G / (ħc * 1000.0)^2)
@printf("  K = %.6f fm⁵ (%.2f MeV⁻⁵)\n", K, K / (ħc * 1000.0)^5)
println()

println("=" ^ 80)
println("Julia 的解 (需要求解能隙方程)")
println("=" ^ 80)
println()

println("注意: 需要实际运行 Julia 的能隙方程求解器来获取解")
println("这里只是展示如何对比")
println()

println("=" ^ 80)
println("关键问题")
println("=" ^ 80)
println()
println("1. Julia 和 Fortran 使用相同的初始猜测吗?")
println("2. Julia 和 Fortran 使用相同的求解器参数吗?")
println("3. Julia 和 Fortran 的能隙方程形式完全一致吗?")
println("4. 在 T=300 MeV, μ=2 MeV 下,两者的解是否收敛到同一点?")
println()

println("=" ^ 80)
println("下一步行动")
println("=" ^ 80)
println()
println("1. 运行 Julia 的能隙方程求解器,获取在 T=300 MeV, μ=2 MeV 下的解")
println("2. 对比 Julia 和 Fortran 的质量值")
println("3. 如果质量不一致,检查:")
println("   - 能隙方程的形式")
println("   - 初始猜测")
println("   - 求解器参数")
println("   - A 函数的计算")
println()
