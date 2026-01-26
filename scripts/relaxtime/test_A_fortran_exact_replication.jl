#!/usr/bin/env julia
"""
完全复制 Fortran 的 A 函数计算逻辑

Fortran 的实现:
1. 常数项: 使用 y(i), w(i) 数值积分 [0, Λ=3.05], 128 节点
2. 分布函数项: 使用 y1(i), w1(i) 数值积分 [0, 15.0], 128 节点
3. 最后乘以 4

目标: 完全复制 Fortran 的计算,看能否得到相同的结果
"""

using Printf

include("../../src/QuarkDistribution.jl")
include("../../src/integration/GaussLegendre.jl")

using .PNJLQuarkDistributions
using .GaussLegendre

# 完全复制 Fortran 的 A 函数计算
function A_fortran_exact(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64)
    # 生成积分节点 (与 Fortran 完全相同)
    y, w = gauleg(0.0, 3.05, 128)      # 常数项
    y1, w1 = gauleg(0.0, 15.0, 128)    # 分布函数项
    
    arrA = 0.0
    
    # 第一个循环: 常数项 (使用 y, w)
    for i in 1:128
        p = y[i]
        E = sqrt(p^2 + m^2)
        arrA += w[i] * p^2 * (-1.0) / E
    end
    
    # 第二个循环: 分布函数项 (使用 y1, w1)
    for i in 1:128
        p = y1[i]
        E = sqrt(p^2 + m^2)
        # Fortran: fphi(E - μ, T, Φ, Φbar) + fphibar(E + μ, T, Φ, Φbar)
        f_quark = quark_distribution(E, μ, T, Φ, Φbar)
        f_antiquark = antiquark_distribution(E, μ, T, Φ, Φbar)
        arrA += w1[i] * p^2 * (f_quark + f_antiquark) / E
    end
    
    # 最后乘以 4
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

# 从 Fortran 结果读取的质量
m_u = 0.00507  # fm⁻¹
m_s = 0.09493  # fm⁻¹

println("=" ^ 80)
println("完全复制 Fortran 的 A 函数计算")
println("=" ^ 80)
println()
println("测试参数:")
@printf("  T = %.6f fm⁻¹ (%.1f MeV)\n", T, T_MeV)
@printf("  μ = %.6f fm⁻¹ (%.1f MeV)\n", μ, μ_MeV)
@printf("  Φ = %.8f\n", Φ)
@printf("  Φbar = %.8f\n", Φbar)
@printf("  m_u = %.5f fm⁻¹\n", m_u)
@printf("  m_s = %.5f fm⁻¹\n", m_s)
println()

println("=" ^ 80)
println("计算结果")
println("=" ^ 80)
println()

A_u_julia = A_fortran_exact(m_u, μ, T, Φ, Φbar)
A_s_julia = A_fortran_exact(m_s, μ, T, Φ, Φbar)

@printf("Julia 计算 (完全复制 Fortran 逻辑):\n")
@printf("  A_u = %.6f fm⁻²\n", A_u_julia)
@printf("  A_s = %.6f fm⁻²\n", A_s_julia)
println()

@printf("Fortran 结果:\n")
@printf("  A_u = -4.972 fm⁻²\n")
@printf("  A_s = (未知)\n")
println()

@printf("差异:\n")
@printf("  ΔA_u = %.6f fm⁻² (%.2f%%)\n", -4.972 - A_u_julia,
        100.0 * abs(-4.972 - A_u_julia) / 4.972)
println()

println("=" ^ 80)
println("详细分析")
println("=" ^ 80)
println()

# 分别计算常数项和分布函数项
function compute_terms(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64)
    y, w = gauleg(0.0, 3.05, 128)
    y1, w1 = gauleg(0.0, 15.0, 128)
    
    const_term = 0.0
    for i in 1:128
        p = y[i]
        E = sqrt(p^2 + m^2)
        const_term += w[i] * p^2 * (-1.0) / E
    end
    
    dist_term = 0.0
    for i in 1:128
        p = y1[i]
        E = sqrt(p^2 + m^2)
        f_quark = quark_distribution(E, μ, T, Φ, Φbar)
        f_antiquark = antiquark_distribution(E, μ, T, Φ, Φbar)
        dist_term += w1[i] * p^2 * (f_quark + f_antiquark) / E
    end
    
    return const_term, dist_term
end

const_term_u, dist_term_u = compute_terms(m_u, μ, T, Φ, Φbar)

@printf("m_u 的详细分解:\n")
@printf("  常数项: %.10f fm⁻²\n", const_term_u)
@printf("  分布函数项: %.10f fm⁻²\n", dist_term_u)
@printf("  总和 (×4): %.10f fm⁻²\n", 4.0 * (const_term_u + dist_term_u))
println()

@printf("  常数项 (×4): %.6f fm⁻²\n", 4.0 * const_term_u)
@printf("  分布函数项 (×4): %.6f fm⁻²\n", 4.0 * dist_term_u)
@printf("  A_u = %.6f fm⁻²\n", A_u_julia)
println()

# 检查分布函数的值
println("=" ^ 80)
println("分布函数值检查 (在几个典型动量点)")
println("=" ^ 80)
println()

test_p_values = [0.5, 1.0, 2.0, 5.0, 10.0]
println("p (fm⁻¹) | E (fm⁻¹) | f_quark | f_antiquark | f_q + f_aq")
println("-" ^ 80)

for p in test_p_values
    E = sqrt(p^2 + m_u^2)
    f_q = quark_distribution(E, μ, T, Φ, Φbar)
    f_aq = antiquark_distribution(E, μ, T, Φ, Φbar)
    @printf("%8.4f | %8.4f | %.6e | %.6e | %.6e\n", p, E, f_q, f_aq, f_q + f_aq)
end

println()
println("=" ^ 80)
println("结论")
println("=" ^ 80)
println()
println("如果 Julia 完全复制 Fortran 的逻辑后仍有差异,")
println("则问题可能在于:")
println("  1. 分布函数的实现细节")
println("  2. 数值精度问题")
println("  3. Fortran 代码中有其他隐藏的处理")
println()
