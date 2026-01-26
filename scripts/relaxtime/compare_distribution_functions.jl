#!/usr/bin/env julia
"""
对比 Fortran 和 Julia 的 PNJL 分布函数

直接对比相同输入下的分布函数值，找出 31% A 函数差异的根源
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "QuarkDistribution.jl"))

using .Constants_PNJL
using .PNJLQuarkDistributions
using Printf

println("="^80)
println("对比 Fortran 和 Julia 的 PNJL 分布函数")
println("="^80)
println()

# 测试参数
T_MeV = 300.0
μ_MeV = 2.0
T = T_MeV / ħc_MeV_fm
μ = μ_MeV / ħc_MeV_fm

Φ = 0.99999994
Φbar = 0.99999994

@printf("物理参数:\n")
@printf("  T = %.1f MeV = %.6f fm⁻¹\n", T_MeV, T)
@printf("  μ = %.1f MeV = %.6f fm⁻¹\n", μ_MeV, μ)
@printf("  Φ = %.8f\n", Φ)
@printf("  Φ̄ = %.8f\n", Φbar)
println()

# Fortran 分布函数的 Julia 实现（根据 Fortran 代码）
function fphi_fortran(x::Float64, T::Float64, Phi1::Float64, Phi2::Float64)
    ee = exp(-x/T)
    return (Phi1*ee + 2.0*Phi2*ee^2 + ee^3) / 
           (1.0 + 3.0*Phi1*ee + 3.0*Phi2*ee^2 + ee^3)
end

function fphibar_fortran(x::Float64, T::Float64, Phi1::Float64, Phi2::Float64)
    ee = exp(-x/T)
    return (Phi2*ee + 2.0*Phi1*ee^2 + ee^3) / 
           (1.0 + 3.0*Phi2*ee + 3.0*Phi1*ee^2 + ee^3)
end

# 测试不同的能量值
println("="^80)
println("测试 1: 不同能量下的分布函数值")
println("="^80)
println()

# 使用 Fortran 的质量
m_u = 0.040510  # fm⁻¹
m_s = 1.032296  # fm⁻¹

# 测试几个代表性的动量值
test_momenta = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0]  # fm⁻¹

println("u 夸克 (m = $(m_u) fm⁻¹):")
println("-"^80)
@printf("%-10s %-15s %-20s %-20s %-15s\n", "p (fm⁻¹)", "E (fm⁻¹)", "f⁺ (Julia)", "f⁺ (Fortran)", "差异 (%)")
println("-"^80)

for p in test_momenta
    E = sqrt(p^2 + m_u^2)
    
    # Julia 的分布函数
    f_plus_julia = quark_distribution(E, μ, T, Φ, Φbar)
    f_minus_julia = antiquark_distribution(E, μ, T, Φ, Φbar)
    
    # Fortran 的分布函数（手动计算）
    x_plus = E - μ
    x_minus = E + μ
    f_plus_fortran = fphi_fortran(x_plus, T, Φ, Φbar)
    f_minus_fortran = fphibar_fortran(x_minus, T, Φ, Φbar)
    
    # 对比 f⁺
    diff_plus = abs(f_plus_julia - f_plus_fortran) / f_plus_fortran * 100
    @printf("%-10.2f %-15.6f %-20.10f %-20.10f %-15.4f\n", 
        p, E, f_plus_julia, f_plus_fortran, diff_plus)
end
println()

println("反夸克分布函数 f⁻:")
println("-"^80)
@printf("%-10s %-15s %-20s %-20s %-15s\n", "p (fm⁻¹)", "E (fm⁻¹)", "f⁻ (Julia)", "f⁻ (Fortran)", "差异 (%)")
println("-"^80)

for p in test_momenta
    E = sqrt(p^2 + m_u^2)
    
    # Julia 的分布函数
    f_minus_julia = antiquark_distribution(E, μ, T, Φ, Φbar)
    
    # Fortran 的分布函数
    x_minus = E + μ
    f_minus_fortran = fphibar_fortran(x_minus, T, Φ, Φbar)
    
    # 对比
    diff_minus = abs(f_minus_julia - f_minus_fortran) / f_minus_fortran * 100
    @printf("%-10.2f %-15.6f %-20.10f %-20.10f %-15.4f\n", 
        p, E, f_minus_julia, f_minus_fortran, diff_minus)
end
println()

# 测试 2: 计算分布函数项的积分贡献
println("="^80)
println("测试 2: 分布函数项对 A 函数的贡献")
println("="^80)
println()

# 使用简单的数值积分估算
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
using .GaussLegendre

nodes_p, weights_p = gauleg(0.0, 15.0, 128)

# Julia 的分布函数积分
integral_julia = 0.0
for i in eachindex(nodes_p)
    p = nodes_p[i]
    w = weights_p[i]
    E = sqrt(p^2 + m_u^2)
    f_plus = quark_distribution(E, μ, T, Φ, Φbar)
    f_minus = antiquark_distribution(E, μ, T, Φ, Φbar)
    integral_julia += w * p^2 / E * (f_plus + f_minus)
end

# Fortran 的分布函数积分
integral_fortran = 0.0
for i in eachindex(nodes_p)
    p = nodes_p[i]
    w = weights_p[i]
    E = sqrt(p^2 + m_u^2)
    x_plus = E - μ
    x_minus = E + μ
    f_plus = fphi_fortran(x_plus, T, Φ, Φbar)
    f_minus = fphibar_fortran(x_minus, T, Φ, Φbar)
    integral_fortran += w * p^2 / E * (f_plus + f_minus)
end

println("分布函数项积分 (m_u, [0, 15 fm⁻¹]):")
@printf("  Julia:   %.6f fm⁻²\n", integral_julia)
@printf("  Fortran: %.6f fm⁻²\n", integral_fortran)
@printf("  差异:    %.2f%%\n", abs(integral_julia - integral_fortran) / integral_fortran * 100)
println()

# 计算完整的 A 值
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
using .OneLoopIntegrals

const_term = OneLoopIntegrals.const_integral_term_A(m_u)

A_julia_manual = 4.0 * (-const_term + integral_julia)
A_fortran_manual = 4.0 * (-const_term + integral_fortran)

println("完整的 A 值 (手动计算):")
@printf("  常数项:  %.6f fm⁻²\n", const_term)
@printf("  A (Julia):   %.6f fm⁻²\n", A_julia_manual)
@printf("  A (Fortran): %.6f fm⁻²\n", A_fortran_manual)
@printf("  差异:        %.2f%%\n", abs(A_julia_manual - A_fortran_manual) / abs(A_fortran_manual) * 100)
println()

# 对比实际的 A 函数结果
A_julia_actual = OneLoopIntegrals.A(m_u, μ, T, Φ, Φbar, nodes_p, weights_p)
A_fortran_actual = -4.972144

println("实际的 A 值:")
@printf("  A (Julia 实际):   %.6f fm⁻²\n", A_julia_actual)
@printf("  A (Fortran 实际): %.6f fm⁻²\n", A_fortran_actual)
@printf("  差异:             %.2f%%\n", abs(A_julia_actual - A_fortran_actual) / abs(A_fortran_actual) * 100)
println()

# 测试 3: 检查分布函数公式
println("="^80)
println("测试 3: 验证分布函数公式")
println("="^80)
println()

# 测试一个具体的点
E_test = 1.0  # fm⁻¹
x_test = E_test - μ

println("测试点: E = $(E_test) fm⁻¹, x = E - μ = $(x_test) fm⁻¹")
println()

# 手动计算 Fortran 公式
ee = exp(-x_test/T)
numerator = Φ*ee + 2.0*Φbar*ee^2 + ee^3
denominator = 1.0 + 3.0*Φ*ee + 3.0*Φbar*ee^2 + ee^3
f_manual = numerator / denominator

println("Fortran 公式 fphi(x, T, Φ, Φ̄):")
@printf("  ee = exp(-x/T) = %.10f\n", ee)
@printf("  分子 = Φ*ee + 2*Φ̄*ee² + ee³ = %.10f\n", numerator)
@printf("  分母 = 1 + 3*Φ*ee + 3*Φ̄*ee² + ee³ = %.10f\n", denominator)
@printf("  f = 分子/分母 = %.10f\n", f_manual)
println()

# Julia 的结果
f_julia = quark_distribution(E_test, μ, T, Φ, Φbar)
println("Julia 的结果:")
@printf("  f⁺(E, μ, T, Φ, Φ̄) = %.10f\n", f_julia)
println()

# 对比
println("对比:")
@printf("  差异 = %.2e (%.4f%%)\n", abs(f_julia - f_manual), abs(f_julia - f_manual) / f_manual * 100)
println()

# 总结
println("="^80)
println("总结")
println("="^80)
println()

max_diff = 0.0
for p in test_momenta
    E = sqrt(p^2 + m_u^2)
    f_plus_julia = quark_distribution(E, μ, T, Φ, Φbar)
    x_plus = E - μ
    f_plus_fortran = fphi_fortran(x_plus, T, Φ, Φbar)
    diff = abs(f_plus_julia - f_plus_fortran) / f_plus_fortran * 100
    max_diff = max(max_diff, diff)
end

if max_diff < 0.1
    println("✅ Julia 和 Fortran 的分布函数完全一致（差异 < 0.1%）")
    println("   → 分布函数不是 31% 差异的原因")
    println()
    println("   需要进一步检查：")
    println("   1. 积分权重的使用")
    println("   2. 能量 E 的计算")
    println("   3. 是否有其他隐藏的实现差异")
elseif max_diff < 1.0
    println("✅ Julia 和 Fortran 的分布函数基本一致（差异 < 1%）")
    println("   → 分布函数可能不是主要原因")
    println("   → 但累积效应可能导致 A 函数差异")
else
    println("❌ Julia 和 Fortran 的分布函数有显著差异（> 1%）")
    println("   → 这可能是 31% A 函数差异的主要原因")
    println()
    println("   需要修正 Julia 的分布函数实现")
end
println()

# 检查积分贡献的差异
integral_diff = abs(integral_julia - integral_fortran) / integral_fortran * 100
if integral_diff > 10.0
    println("⚠️  分布函数项积分差异很大（> 10%）")
    println("   → 这可以解释 31% 的 A 函数差异")
elseif integral_diff > 1.0
    println("⚠️  分布函数项积分有一定差异（> 1%）")
    println("   → 可能部分解释 A 函数差异")
else
    println("✅ 分布函数项积分基本一致（< 1%）")
    println("   → 需要寻找其他原因")
end
println()

println("="^80)
