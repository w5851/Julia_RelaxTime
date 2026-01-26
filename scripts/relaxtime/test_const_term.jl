#!/usr/bin/env julia
"""
测试常数项的计算

对比 Julia 的解析公式和 Fortran 的数值积分
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
println("测试常数项的计算")
println("="^80)
println()

# 测试参数
m_u = 0.040510  # fm⁻¹
m_s = 1.032296  # fm⁻¹

println("测试质量:")
@printf("  m_u = %.6f fm⁻¹\n", m_u)
@printf("  m_s = %.6f fm⁻¹\n", m_s)
println()

# Julia 的解析公式
const_term_u_julia = OneLoopIntegrals.const_integral_term_A(m_u)
const_term_s_julia = OneLoopIntegrals.const_integral_term_A(m_s)

println("Julia 解析公式:")
@printf("  ∫₀^Λ p²/E dp (m_u) = %.6f fm⁻²\n", const_term_u_julia)
@printf("  ∫₀^Λ p²/E dp (m_s) = %.6f fm⁻²\n", const_term_s_julia)
println()

# Fortran 的数值积分（模拟）
Λ = Λ_inv_fm  # 3.0523 fm⁻¹
nodes_Lambda, weights_Lambda = gauleg(0.0, Λ, 128)

const_term_u_numerical = sum(weights_Lambda[i] * nodes_Lambda[i]^2 / sqrt(nodes_Lambda[i]^2 + m_u^2) 
                             for i in eachindex(nodes_Lambda))
const_term_s_numerical = sum(weights_Lambda[i] * nodes_Lambda[i]^2 / sqrt(nodes_Lambda[i]^2 + m_s^2) 
                             for i in eachindex(nodes_Lambda))

println("Fortran 数值积分（模拟，128 节点）:")
@printf("  ∫₀^Λ p²/E dp (m_u) = %.6f fm⁻²\n", const_term_u_numerical)
@printf("  ∫₀^Λ p²/E dp (m_s) = %.6f fm⁻²\n", const_term_s_numerical)
println()

# 对比
println("差异:")
diff_u = abs(const_term_u_julia - const_term_u_numerical) / const_term_u_numerical * 100
diff_s = abs(const_term_s_julia - const_term_s_numerical) / const_term_s_numerical * 100
@printf("  Δ(m_u) = %.4f%% (%.6e fm⁻²)\n", diff_u, abs(const_term_u_julia - const_term_u_numerical))
@printf("  Δ(m_s) = %.4f%% (%.6e fm⁻²)\n", diff_s, abs(const_term_s_julia - const_term_s_numerical))
println()

# 检查公式
println("="^80)
println("验证解析公式")
println("="^80)
println()

println("解析公式：∫₀^Λ p²/E dp = (Λ/2)√(Λ²+m²) - (m²/2)ln((Λ+√(Λ²+m²))/m)")
println()

# 手动计算
sqrt_term_u = sqrt(Λ^2 + m_u^2)
term1_u = (Λ / 2.0) * sqrt_term_u
term2_u = (m_u^2 / 2.0) * log((Λ + sqrt_term_u) / m_u)
manual_u = term1_u - term2_u

sqrt_term_s = sqrt(Λ^2 + m_s^2)
term1_s = (Λ / 2.0) * sqrt_term_s
term2_s = (m_s^2 / 2.0) * log((Λ + sqrt_term_s) / m_s)
manual_s = term1_s - term2_s

println("手动计算 (m_u):")
@printf("  term1 = (Λ/2)√(Λ²+m²) = %.6f\n", term1_u)
@printf("  term2 = (m²/2)ln(...) = %.6f\n", term2_u)
@printf("  结果 = term1 - term2 = %.6f fm⁻²\n", manual_u)
println()

println("手动计算 (m_s):")
@printf("  term1 = (Λ/2)√(Λ²+m²) = %.6f\n", term1_s)
@printf("  term2 = (m²/2)ln(...) = %.6f\n", term2_s)
@printf("  结果 = term1 - term2 = %.6f fm⁻²\n", manual_s)
println()

# 总结
println("="^80)
println("总结")
println("="^80)
println()

if diff_u < 0.1 && diff_s < 0.1
    println("✅ Julia 的解析公式与数值积分一致（差异 < 0.1%）")
    println("   → 常数项的计算是正确的")
    println("   → 31% 的差异不是来自常数项")
else
    println("❌ Julia 的解析公式与数值积分有差异（> 0.1%）")
    println("   → 需要检查解析公式的实现")
end
println()

# 估算常数项对总 A 值的贡献
println("="^80)
println("常数项对总 A 值的贡献")
println("="^80)
println()

# Fortran 的 A 值
A_u_fortran = -4.972144
A_s_fortran = -5.165541

# 常数项的贡献（乘以 -4）
const_contrib_u = -4.0 * const_term_u_julia
const_contrib_s = -4.0 * const_term_s_julia

println("常数项贡献（乘以 -4）:")
@printf("  m_u: %.6f fm⁻² (占 Fortran A_u 的 %.1f%%)\n", 
    const_contrib_u, abs(const_contrib_u / A_u_fortran) * 100)
@printf("  m_s: %.6f fm⁻² (占 Fortran A_s 的 %.1f%%)\n", 
    const_contrib_s, abs(const_contrib_s / A_s_fortran) * 100)
println()

println("="^80)
