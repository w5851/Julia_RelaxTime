#!/usr/bin/env julia
"""
验证 Polyakov 环假设

使用 Fortran 的 Φ 值重新计算 Julia 的 A 函数,
看差异是否消失
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
const ħc = 0.1973269804
T_MeV = 300.0
μ_MeV = 2.0
T = T_MeV / (ħc * 1000.0)
μ = μ_MeV / (ħc * 1000.0)

# Fortran 的参数
m_u_fortran = 0.04051  # fm⁻¹
m_s_fortran = 1.0323   # fm⁻¹
Phi_fortran = 0.99999994
Phibar_fortran = 0.99999994

# Julia 的参数
m_u_julia = 0.04122  # fm⁻¹
m_s_julia = 1.04815  # fm⁻¹
Phi_julia = 0.838
Phibar_julia = 0.838

println("=" ^ 80)
println("验证 Polyakov 环假设")
println("=" ^ 80)
println()

println("=" ^ 80)
println("测试 1: 使用 Julia 的质量 + Julia 的 Φ")
println("=" ^ 80)
println()

A_u_julia_julia = A_fortran_exact(m_u_julia, μ, T, Phi_julia, Phibar_julia)
A_s_julia_julia = A_fortran_exact(m_s_julia, μ, T, Phi_julia, Phibar_julia)

@printf("m_u = %.5f fm⁻¹, Φ = %.3f:\n", m_u_julia, Phi_julia)
@printf("  A_u = %.6f fm⁻²\n", A_u_julia_julia)
println()

@printf("m_s = %.5f fm⁻¹, Φ = %.3f:\n", m_s_julia, Phi_julia)
@printf("  A_s = %.6f fm⁻²\n", A_s_julia_julia)
println()

println("=" ^ 80)
println("测试 2: 使用 Fortran 的质量 + Fortran 的 Φ")
println("=" ^ 80)
println()

A_u_fortran_fortran = A_fortran_exact(m_u_fortran, μ, T, Phi_fortran, Phibar_fortran)
A_s_fortran_fortran = A_fortran_exact(m_s_fortran, μ, T, Phi_fortran, Phibar_fortran)

@printf("m_u = %.5f fm⁻¹, Φ = %.3f:\n", m_u_fortran, Phi_fortran)
@printf("  A_u = %.6f fm⁻²\n", A_u_fortran_fortran)
println()

@printf("m_s = %.5f fm⁻¹, Φ = %.3f:\n", m_s_fortran, Phi_fortran)
@printf("  A_s = %.6f fm⁻²\n", A_s_fortran_fortran)
println()

println("=" ^ 80)
println("测试 3: 使用 Julia 的质量 + Fortran 的 Φ (关键测试!)")
println("=" ^ 80)
println()

A_u_julia_fortran = A_fortran_exact(m_u_julia, μ, T, Phi_fortran, Phibar_fortran)
A_s_julia_fortran = A_fortran_exact(m_s_julia, μ, T, Phi_fortran, Phibar_fortran)

@printf("m_u = %.5f fm⁻¹, Φ = %.3f:\n", m_u_julia, Phi_fortran)
@printf("  A_u = %.6f fm⁻²\n", A_u_julia_fortran)
println()

@printf("m_s = %.5f fm⁻¹, Φ = %.3f:\n", m_s_julia, Phi_fortran)
@printf("  A_s = %.6f fm⁻²\n", A_s_julia_fortran)
println()

println("=" ^ 80)
println("测试 4: 使用 Fortran 的质量 + Julia 的 Φ")
println("=" ^ 80)
println()

A_u_fortran_julia = A_fortran_exact(m_u_fortran, μ, T, Phi_julia, Phibar_julia)
A_s_fortran_julia = A_fortran_exact(m_s_fortran, μ, T, Phi_julia, Phibar_julia)

@printf("m_u = %.5f fm⁻¹, Φ = %.3f:\n", m_u_fortran, Phi_julia)
@printf("  A_u = %.6f fm⁻²\n", A_u_fortran_julia)
println()

@printf("m_s = %.5f fm⁻¹, Φ = %.3f:\n", m_s_fortran, Phi_julia)
@printf("  A_s = %.6f fm⁻²\n", A_s_fortran_julia)
println()

println("=" ^ 80)
println("对比 Fortran 结果")
println("=" ^ 80)
println()

A_u_fortran_target = -4.972162
A_s_fortran_target = -5.165541

@printf("Fortran 目标值:\n")
@printf("  A_u = %.6f fm⁻²\n", A_u_fortran_target)
@printf("  A_s = %.6f fm⁻²\n", A_s_fortran_target)
println()

println("=" ^ 80)
println("差异分析")
println("=" ^ 80)
println()

@printf("测试 2 (Fortran 质量 + Fortran Φ) vs Fortran 目标:\n")
@printf("  ΔA_u = %.6f fm⁻² (%.2f%%)\n", 
        A_u_fortran_target - A_u_fortran_fortran,
        100.0 * abs(A_u_fortran_target - A_u_fortran_fortran) / abs(A_u_fortran_target))
@printf("  ΔA_s = %.6f fm⁻² (%.2f%%)\n",
        A_s_fortran_target - A_s_fortran_fortran,
        100.0 * abs(A_s_fortran_target - A_s_fortran_fortran) / abs(A_s_fortran_target))
println()

@printf("测试 3 (Julia 质量 + Fortran Φ) vs Fortran 目标:\n")
@printf("  ΔA_u = %.6f fm⁻² (%.2f%%)\n",
        A_u_fortran_target - A_u_julia_fortran,
        100.0 * abs(A_u_fortran_target - A_u_julia_fortran) / abs(A_u_fortran_target))
@printf("  ΔA_s = %.6f fm⁻² (%.2f%%)\n",
        A_s_fortran_target - A_s_julia_fortran,
        100.0 * abs(A_s_fortran_target - A_s_julia_fortran) / abs(A_s_fortran_target))
println()

println("=" ^ 80)
println("Φ 的影响分析")
println("=" ^ 80)
println()

@printf("固定 Fortran 质量,改变 Φ:\n")
@printf("  Φ = %.3f: A_u = %.6f fm⁻²\n", Phi_julia, A_u_fortran_julia)
@printf("  Φ = %.3f: A_u = %.6f fm⁻²\n", Phi_fortran, A_u_fortran_fortran)
@printf("  差异: %.6f fm⁻² (%.2f%%)\n",
        A_u_fortran_fortran - A_u_fortran_julia,
        100.0 * abs(A_u_fortran_fortran - A_u_fortran_julia) / abs(A_u_fortran_julia))
println()

println("=" ^ 80)
println("结论")
println("=" ^ 80)
println()

if abs(A_u_fortran_target - A_u_fortran_fortran) / abs(A_u_fortran_target) < 0.02
    println("✅ 使用 Fortran 的质量和 Φ,Julia 可以复现 Fortran 的 A 函数!")
    println()
    println("这证明:")
    println("  1. Julia 的 A 函数实现是正确的")
    println("  2. 31% 的差异完全来自 Polyakov 环 Φ 的不同")
    println("  3. Julia 和 Fortran 求解到了不同的能隙方程解")
    println()
    println("根本原因: Polyakov 环的差异 (Julia: 0.838 vs Fortran: 1.0)")
else
    println("❌ 即使使用 Fortran 的参数,仍有差异")
    println()
    println("需要进一步检查:")
    println("  1. A 函数的实现细节")
    println("  2. 分布函数的实现")
    println("  3. 其他隐藏的差异")
end
println()
