#!/usr/bin/env julia
"""
详细对比 Julia 和 Fortran 的 A 函数计算

测试点: T=300 MeV, μ=2 MeV
- Fortran A_u = -4.972 fm⁻²
- Julia A_u = -3.426 fm⁻² (使用 20 fm⁻¹ 积分上限)
- 差异: 31%

目标: 找出差异来源
"""

using Printf

# 包含必要的模块
include("../../src/QuarkDistribution.jl")
include("../../src/Constants_PNJL.jl")
include("../../src/integration/GaussLegendre.jl")

using .PNJLQuarkDistributions
using .Constants_PNJL
using .GaussLegendre

# 常数项积分的解析公式
function const_integral_term_A(m::Float64, Lambda::Float64)
    m_pos = max(m, 0.0)
    if m_pos < 1e-14
        return (Lambda^2) / 2.0
    end
    term1 = Lambda * sqrt(Lambda^2 + m_pos^2)
    term2 = m_pos^2 * log((Lambda + sqrt(Lambda^2 + m_pos^2)) / m_pos)
    return (term1 - term2) / 2.0
end

# Fortran 风格的 A 函数计算
function A_fortran_style(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64,
                         Lambda_const::Float64, Lambda_dist::Float64, n_points::Int)
    # 常数项: 使用 Lambda_const 截断
    const_term = -const_integral_term_A(m, Lambda_const)
    
    # 分布函数项: 使用 Lambda_dist 积分
    nodes, weights = gauleg(0.0, Lambda_dist, n_points)
    dist_term = 0.0
    for i in 1:n_points
        p = nodes[i]
        w = weights[i]
        E = sqrt(p^2 + m^2)
        f_quark = quark_distribution(E, μ, T, Φ, Φbar)
        f_antiquark = antiquark_distribution(E, μ, T, Φ, Φbar)
        dist_term += w * p^2 / E * (f_quark + f_antiquark)
    end
    
    return 4.0 * (const_term + dist_term), const_term, dist_term
end

# Julia 风格的 A 函数计算(当前实现)
function A_julia_style(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64,
                       Lambda_const::Float64, Lambda_dist::Float64, n_points::Int)
    # 与 Fortran 风格相同,只是为了对比
    return A_fortran_style(m, μ, T, Φ, Φbar, Lambda_const, Lambda_dist, n_points)
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
println("A 函数详细对比测试")
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

# Fortran 使用的积分范围
Lambda_const_fortran = 3.05  # fm⁻¹ (Λ截断)
Lambda_dist_fortran = 15.0   # fm⁻¹ (分布函数积分上限)

# Julia 当前使用的积分范围
Lambda_const_julia = 3.05    # fm⁻¹ (Λ截断)
Lambda_dist_julia = 20.0     # fm⁻¹ (分布函数积分上限)

println("=" ^ 80)
println("1. 常数项对比")
println("=" ^ 80)
println()

const_u_fortran = const_integral_term_A(m_u, Lambda_const_fortran)
const_s_fortran = const_integral_term_A(m_s, Lambda_const_fortran)
const_u_julia = const_integral_term_A(m_u, Lambda_const_julia)
const_s_julia = const_integral_term_A(m_s, Lambda_const_julia)

@printf("m_u 常数项:\n")
@printf("  Fortran (Λ=%.2f): %.10f fm⁻²\n", Lambda_const_fortran, const_u_fortran)
@printf("  Julia   (Λ=%.2f): %.10f fm⁻²\n", Lambda_const_julia, const_u_julia)
@printf("  差异: %.6e\n", abs(const_u_fortran - const_u_julia))
println()

@printf("m_s 常数项:\n")
@printf("  Fortran (Λ=%.2f): %.10f fm⁻²\n", Lambda_const_fortran, const_s_fortran)
@printf("  Julia   (Λ=%.2f): %.10f fm⁻²\n", Lambda_const_julia, const_s_julia)
@printf("  差异: %.6e\n", abs(const_s_fortran - const_s_julia))
println()

println("=" ^ 80)
println("2. 完整 A 函数对比 (使用不同积分上限)")
println("=" ^ 80)
println()

# 测试不同的积分节点数
n_points_list = [16, 32, 64, 128]

for n_points in n_points_list
    println("-" ^ 80)
    @printf("使用 %d 个积分节点:\n", n_points)
    println()
    
    # Fortran 风格 (Λ_const=3.05, Λ_dist=15.0)
    A_u_f, const_u_f, dist_u_f = A_fortran_style(m_u, μ, T, Φ, Φbar, 
                                                   Lambda_const_fortran, Lambda_dist_fortran, n_points)
    A_s_f, const_s_f, dist_s_f = A_fortran_style(m_s, μ, T, Φ, Φbar,
                                                   Lambda_const_fortran, Lambda_dist_fortran, n_points)
    
    # Julia 风格 (Λ_const=3.05, Λ_dist=20.0)
    A_u_j, const_u_j, dist_u_j = A_julia_style(m_u, μ, T, Φ, Φbar,
                                                 Lambda_const_julia, Lambda_dist_julia, n_points)
    A_s_j, const_s_j, dist_s_j = A_julia_style(m_s, μ, T, Φ, Φbar,
                                                 Lambda_const_julia, Lambda_dist_julia, n_points)
    
    @printf("A_u:\n")
    @printf("  Fortran (15 fm⁻¹): %.6f fm⁻² (常数: %.6f, 分布: %.6f)\n", 
            A_u_f, 4.0*const_u_f, 4.0*dist_u_f)
    @printf("  Julia   (20 fm⁻¹): %.6f fm⁻² (常数: %.6f, 分布: %.6f)\n",
            A_u_j, 4.0*const_u_j, 4.0*dist_u_j)
    @printf("  差异: %.6f fm⁻² (%.2f%%)\n", A_u_f - A_u_j, 
            100.0 * abs(A_u_f - A_u_j) / abs(A_u_f))
    println()
    
    @printf("A_s:\n")
    @printf("  Fortran (15 fm⁻¹): %.6f fm⁻² (常数: %.6f, 分布: %.6f)\n",
            A_s_f, 4.0*const_s_f, 4.0*dist_s_f)
    @printf("  Julia   (20 fm⁻¹): %.6f fm⁻² (常数: %.6f, 分布: %.6f)\n",
            A_s_j, 4.0*const_s_j, 4.0*dist_s_j)
    @printf("  差异: %.6f fm⁻² (%.2f%%)\n", A_s_f - A_s_j,
            100.0 * abs(A_s_f - A_s_j) / abs(A_s_f))
    println()
end

println("=" ^ 80)
println("3. 测试使用相同积分上限 (15 fm⁻¹)")
println("=" ^ 80)
println()

n_points = 64
A_u_same, const_u_same, dist_u_same = A_fortran_style(m_u, μ, T, Φ, Φbar,
                                                        Lambda_const_fortran, Lambda_dist_fortran, n_points)
A_s_same, const_s_same, dist_s_same = A_fortran_style(m_s, μ, T, Φ, Φbar,
                                                        Lambda_const_fortran, Lambda_dist_fortran, n_points)

@printf("使用相同参数 (Λ_const=3.05, Λ_dist=15.0, n=%d):\n", n_points)
@printf("  A_u = %.6f fm⁻²\n", A_u_same)
@printf("  A_s = %.6f fm⁻²\n", A_s_same)
println()
@printf("与 Fortran 结果对比:\n")
@printf("  Fortran A_u = -4.972 fm⁻²\n")
@printf("  Julia   A_u = %.6f fm⁻²\n", A_u_same)
@printf("  差异: %.6f fm⁻² (%.2f%%)\n", -4.972 - A_u_same,
        100.0 * abs(-4.972 - A_u_same) / 4.972)
println()

println("=" ^ 80)
println("结论")
println("=" ^ 80)
println()
println("如果使用相同积分上限后差异消失,则问题在于积分范围")
println("如果差异仍然存在,则需要检查:")
println("  1. 积分节点的生成方式")
println("  2. 被积函数的计算")
println("  3. 其他隐藏的实现差异")
println()
