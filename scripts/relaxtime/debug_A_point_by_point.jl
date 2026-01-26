#!/usr/bin/env julia
"""
逐点对比 A 函数的计算

目标: 找出 Julia 和 Fortran 在哪个环节产生了差异
"""

using Printf

include("../../src/QuarkDistribution.jl")
include("../../src/integration/GaussLegendre.jl")

using .PNJLQuarkDistributions
using .GaussLegendre

# 测试参数
const ħc = 0.1973269804
T_MeV = 300.0
μ_MeV = 2.0
T = T_MeV / (ħc * 1000.0)
μ = μ_MeV / (ħc * 1000.0)
Φ = 0.99999994
Φbar = 0.99999994
m_u = 0.04051  # fm⁻¹

println("=" ^ 80)
println("逐点调试 A 函数计算")
println("=" ^ 80)
println()

# 生成积分节点
y, w = gauleg(0.0, 3.05, 128)
y1, w1 = gauleg(0.0, 15.0, 128)

println("检查前几个积分节点:")
println()
println("常数项积分节点 (y, w):")
println("i | p (fm⁻¹) | weight | E (fm⁻¹) | p²/E | contribution")
println("-" ^ 80)
for i in 1:5
    p = y[i]
    weight = w[i]
    E = sqrt(p^2 + m_u^2)
    integrand = p^2 / E
    contrib = weight * integrand * (-1.0)
    @printf("%d | %.6f | %.6e | %.6f | %.6f | %.6e\n", i, p, weight, E, integrand, contrib)
end
println("...")
println()

println("分布函数项积分节点 (y1, w1):")
println("i | p (fm⁻¹) | weight | E (fm⁻¹) | f_q | f_aq | f_q+f_aq | p²/E×(f_q+f_aq) | contribution")
println("-" ^ 100)
for i in 1:5
    p = y1[i]
    weight = w1[i]
    E = sqrt(p^2 + m_u^2)
    f_q = quark_distribution(E, μ, T, Φ, Φbar)
    f_aq = antiquark_distribution(E, μ, T, Φ, Φbar)
    f_sum = f_q + f_aq
    integrand = p^2 / E * f_sum
    contrib = weight * integrand
    @printf("%d | %.6f | %.6e | %.6f | %.6e | %.6e | %.6e | %.6e | %.6e\n", 
            i, p, weight, E, f_q, f_aq, f_sum, integrand, contrib)
end
println("...")
println()

# 完整计算
function compute_A(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64)
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
        f_q = quark_distribution(E, μ, T, Φ, Φbar)
        f_aq = antiquark_distribution(E, μ, T, Φ, Φbar)
        dist_term += w1[i] * p^2 * (f_q + f_aq) / E
    end
    
    return const_term, dist_term, 4.0 * (const_term + dist_term)
end

const_term, dist_term, A_u = compute_A(m_u, μ, T, Φ, Φbar)

println("=" ^ 80)
println("完整计算结果")
println("=" ^ 80)
println()
@printf("常数项: %.10f fm⁻²\n", const_term)
@printf("分布函数项: %.10f fm⁻²\n", dist_term)
@printf("总和: %.10f fm⁻²\n", const_term + dist_term)
@printf("A_u (×4): %.10f fm⁻²\n", A_u)
println()
@printf("Fortran A_u: -4.972162 fm⁻²\n")
@printf("差异: %.6f fm⁻² (%.2f%%)\n", -4.972162 - A_u, 100.0 * abs(-4.972162 - A_u) / 4.972162)
println()

# 检查分布函数在高动量区域的行为
println("=" ^ 80)
println("检查分布函数在不同动量区域的行为")
println("=" ^ 80)
println()

p_test = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0]
println("p (fm⁻¹) | E (fm⁻¹) | f_quark | f_antiquark | f_q + f_aq")
println("-" ^ 80)
for p in p_test
    E = sqrt(p^2 + m_u^2)
    f_q = quark_distribution(E, μ, T, Φ, Φbar)
    f_aq = antiquark_distribution(E, μ, T, Φ, Φbar)
    @printf("%8.4f | %8.4f | %.6e | %.6e | %.6e\n", p, E, f_q, f_aq, f_q + f_aq)
end
println()

println("=" ^ 80)
println("可能的问题:")
println("=" ^ 80)
println()
println("1. 分布函数的实现可能有细微差异")
println("2. Fortran 可能使用了不同的 Φ 和 Φbar 值")
println("3. Fortran 可能有额外的归一化因子")
println("4. 积分节点的生成可能有差异")
println()
