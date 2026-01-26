#!/usr/bin/env julia
"""
使用正确的质量重新测试 A 函数

从 Fortran 调试输出:
  m_u = 0.04051 fm⁻¹ (不是 0.00507!)
  m_s = 1.0323 fm⁻¹ (不是 0.09493!)
  A_u = -4.972 fm⁻²

这是一个巨大的差异! 让我们用正确的质量重新计算
"""

using Printf

include("../../src/QuarkDistribution.jl")
include("../../src/integration/GaussLegendre.jl")

using .PNJLQuarkDistributions
using .GaussLegendre

# 完全复制 Fortran 的 A 函数计算
function A_fortran_exact(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64)
    y, w = gauleg(0.0, 3.05, 128)
    y1, w1 = gauleg(0.0, 15.0, 128)
    
    arrA = 0.0
    
    # 常数项
    for i in 1:128
        p = y[i]
        E = sqrt(p^2 + m^2)
        arrA += w[i] * p^2 * (-1.0) / E
    end
    
    # 分布函数项
    for i in 1:128
        p = y1[i]
        E = sqrt(p^2 + m^2)
        f_quark = quark_distribution(E, μ, T, Φ, Φbar)
        f_antiquark = antiquark_distribution(E, μ, T, Φ, Φbar)
        arrA += w1[i] * p^2 * (f_quark + f_antiquark) / E
    end
    
    arrA *= 4.0
    return arrA
end

# 测试参数
const ħc = 0.1973269804  # GeV·fm
T_MeV = 300.0
μ_MeV = 2.0
T = T_MeV / (ħc * 1000.0)
μ = μ_MeV / (ħc * 1000.0)
Φ = 0.99999994
Φbar = 0.99999994

# 使用正确的质量 (从 Fortran 调试输出)
m_u_correct = 0.04051  # fm⁻¹
m_s_correct = 1.0323   # fm⁻¹

println("=" ^ 80)
println("使用正确的质量重新测试 A 函数")
println("=" ^ 80)
println()
println("测试参数:")
@printf("  T = %.6f fm⁻¹ (%.1f MeV)\n", T, T_MeV)
@printf("  μ = %.6f fm⁻¹ (%.1f MeV)\n", μ, μ_MeV)
@printf("  Φ = %.8f\n", Φ)
@printf("  Φbar = %.8f\n", Φbar)
println()
@printf("质量 (从 Fortran 调试输出):\n")
@printf("  m_u = %.5f fm⁻¹ (%.1f MeV)\n", m_u_correct, m_u_correct * ħc * 1000.0)
@printf("  m_s = %.5f fm⁻¹ (%.1f MeV)\n", m_s_correct, m_s_correct * ħc * 1000.0)
println()

println("=" ^ 80)
println("计算结果")
println("=" ^ 80)
println()

A_u_julia = A_fortran_exact(m_u_correct, μ, T, Φ, Φbar)
A_s_julia = A_fortran_exact(m_s_correct, μ, T, Φ, Φbar)

@printf("Julia 计算:\n")
@printf("  A_u = %.6f fm⁻²\n", A_u_julia)
@printf("  A_s = %.6f fm⁻²\n", A_s_julia)
println()

@printf("Fortran 结果 (从调试输出):\n")
@printf("  A_u = -4.972162 fm⁻²\n")
@printf("  A_s = -5.165541 fm⁻²\n")
println()

@printf("差异:\n")
@printf("  ΔA_u = %.6f fm⁻² (%.2f%%)\n", -4.972162 - A_u_julia,
        100.0 * abs(-4.972162 - A_u_julia) / 4.972162)
@printf("  ΔA_s = %.6f fm⁻² (%.2f%%)\n", -5.165541 - A_s_julia,
        100.0 * abs(-5.165541 - A_s_julia) / 5.165541)
println()

println("=" ^ 80)
println("结论")
println("=" ^ 80)
println()
println("如果差异消失,则说明:")
println("  1. 之前使用的质量值不正确")
println("  2. Julia 和 Fortran 的 A 函数实现是一致的")
println("  3. 31% 的差异是由于使用了错误的质量值")
println()
println("如果差异仍然存在,则需要进一步调查")
println()
