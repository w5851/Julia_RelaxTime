#!/usr/bin/env julia
"""
对比常数项的数值积分和解析公式

Fortran 使用数值积分:
  ∫₀^Λ p²/E dp  (使用高斯-勒让德积分)

Julia 使用解析公式:
  Λ/2 × √(Λ² + m²) - m²/2 × ln((Λ + √(Λ² + m²))/m)

测试: 这两种方法的差异有多大?
"""

using Printf

include("../../src/integration/GaussLegendre.jl")
using .GaussLegendre

# 解析公式
function const_term_analytical(m::Float64, Lambda::Float64)
    m_pos = max(m, 0.0)
    if m_pos < 1e-14
        return (Lambda^2) / 2.0
    end
    term1 = Lambda * sqrt(Lambda^2 + m_pos^2)
    term2 = m_pos^2 * log((Lambda + sqrt(Lambda^2 + m_pos^2)) / m_pos)
    return (term1 - term2) / 2.0
end

# 数值积分
function const_term_numerical(m::Float64, Lambda::Float64, n_points::Int)
    nodes, weights = gauleg(0.0, Lambda, n_points)
    integral = 0.0
    for i in 1:n_points
        p = nodes[i]
        w = weights[i]
        E = sqrt(p^2 + m^2)
        integral += w * p^2 / E
    end
    return integral
end

# 测试参数
m_u = 0.00507  # fm⁻¹
m_s = 0.09493  # fm⁻¹
Lambda = 3.05  # fm⁻¹

println("=" ^ 80)
println("常数项: 数值积分 vs 解析公式")
println("=" ^ 80)
println()
println("测试参数:")
@printf("  m_u = %.5f fm⁻¹\n", m_u)
@printf("  m_s = %.5f fm⁻¹\n", m_s)
@printf("  Λ = %.2f fm⁻¹\n", Lambda)
println()

# 测试不同的积分节点数
n_points_list = [16, 32, 64, 128, 256]

println("=" ^ 80)
println("m_u 常数项对比")
println("=" ^ 80)
println()
println("节点数 | 数值积分 | 解析公式 | 相对差异")
println("-" ^ 80)

analytical_u = const_term_analytical(m_u, Lambda)

for n_points in n_points_list
    numerical_u = const_term_numerical(m_u, Lambda, n_points)
    rel_diff = abs(numerical_u - analytical_u) / abs(analytical_u)
    @printf("%6d | %.10f | %.10f | %.6e\n", n_points, numerical_u, analytical_u, rel_diff)
end

println()
println("=" ^ 80)
println("m_s 常数项对比")
println("=" ^ 80)
println()
println("节点数 | 数值积分 | 解析公式 | 相对差异")
println("-" ^ 80)

analytical_s = const_term_analytical(m_s, Lambda)

for n_points in n_points_list
    numerical_s = const_term_numerical(m_s, Lambda, n_points)
    rel_diff = abs(numerical_s - analytical_s) / abs(analytical_s)
    @printf("%6d | %.10f | %.10f | %.6e\n", n_points, numerical_s, analytical_s, rel_diff)
end

println()
println("=" ^ 80)
println("结论")
println("=" ^ 80)
println()
println("如果相对差异 < 1e-10，则数值积分和解析公式一致")
println("Fortran 使用 128 个节点的数值积分")
println()

# 计算 Fortran 使用的值
numerical_u_128 = const_term_numerical(m_u, Lambda, 128)
numerical_s_128 = const_term_numerical(m_s, Lambda, 128)

println("Fortran 使用的常数项 (128 节点):")
@printf("  m_u: %.10f fm⁻²\n", numerical_u_128)
@printf("  m_s: %.10f fm⁻²\n", numerical_s_128)
println()
println("Julia 使用的常数项 (解析公式):")
@printf("  m_u: %.10f fm⁻²\n", analytical_u)
@printf("  m_s: %.10f fm⁻²\n", analytical_s)
println()
