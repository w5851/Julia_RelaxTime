#!/usr/bin/env julia
"""
使用 Fortran 的解作为初始值在 Julia 中求解

目的：验证 Julia 的求解器能否收敛到 Fortran 的解

测试条件：T=300 MeV, μ=2 MeV
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "solver", "Solver.jl"))

using .Constants_PNJL
using .Solver
using Printf

println("="^80)
println("使用 Fortran 初始值测试 Julia 求解器")
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

# Fortran 的解（作为初始值）
println("="^80)
println("Fortran 的解（用作初始值）")
println("="^80)

m_u_fortran = 0.040510  # fm⁻¹
m_d_fortran = 0.040510  # fm⁻¹
m_s_fortran = 1.032296  # fm⁻¹
Φ_fortran = 0.99999994
Φbar_fortran = 0.99999994

@printf("  m_u     = %.6f fm⁻¹ = %.2f MeV\n", m_u_fortran, m_u_fortran * ħc_MeV_fm)
@printf("  m_d     = %.6f fm⁻¹ = %.2f MeV\n", m_d_fortran, m_d_fortran * ħc_MeV_fm)
@printf("  m_s     = %.6f fm⁻¹ = %.2f MeV\n", m_s_fortran, m_s_fortran * ħc_MeV_fm)
@printf("  Φ       = %.8f\n", Φ_fortran)
@printf("  Φ̄       = %.8f\n", Φbar_fortran)
println()

# 测试 1：使用 Fortran 的解作为初始值
println("="^80)
println("测试 1：使用 Fortran 解作为初始值")
println("="^80)

# 创建 Fortran 初始值的种子策略
fortran_seed_vec = [m_u_fortran, m_d_fortran, m_s_fortran, Φ_fortran, Φbar_fortran]

# 创建自定义种子策略
using StaticArrays
struct FortranSeed <: SeedStrategy 
    seed::Vector{Float64}
end
Solver.SeedStrategies.get_seed(fs::FortranSeed, ::AbstractVector, ::ConstraintMode) = SVector{5, Float64}(fs.seed...)

result1 = solve(FixedMu(), T, μ; seed_strategy=FortranSeed(fortran_seed_vec))

println("\nJulia 求解结果（Fortran 初始值）:")
@printf("  收敛状态: %s\n", result1.converged ? "✅ 收敛" : "❌ 未收敛")
@printf("  m_u     = %.6f fm⁻¹ = %.2f MeV\n", result1.masses[1], result1.masses[1] * ħc_MeV_fm)
@printf("  m_d     = %.6f fm⁻¹ = %.2f MeV\n", result1.masses[2], result1.masses[2] * ħc_MeV_fm)
@printf("  m_s     = %.6f fm⁻¹ = %.2f MeV\n", result1.masses[3], result1.masses[3] * ħc_MeV_fm)
@printf("  Φ       = %.8f\n", result1.Phi)
@printf("  Φ̄       = %.8f\n", result1.Phibar)
println()

# 对比
println("与 Fortran 解的差异:")
@printf("  Δm_u    = %.2e fm⁻¹ (%.4f%%)\n", 
    abs(result1.masses[1] - m_u_fortran), 
    abs(result1.masses[1] - m_u_fortran) / m_u_fortran * 100)
@printf("  Δm_d    = %.2e fm⁻¹ (%.4f%%)\n", 
    abs(result1.masses[2] - m_d_fortran), 
    abs(result1.masses[2] - m_d_fortran) / m_d_fortran * 100)
@printf("  Δm_s    = %.2e fm⁻¹ (%.4f%%)\n", 
    abs(result1.masses[3] - m_s_fortran), 
    abs(result1.masses[3] - m_s_fortran) / m_s_fortran * 100)
@printf("  ΔΦ      = %.2e (%.4f%%)\n", 
    abs(result1.Phi - Φ_fortran), 
    abs(result1.Phi - Φ_fortran) / Φ_fortran * 100)
@printf("  ΔΦ̄      = %.2e (%.4f%%)\n", 
    abs(result1.Phibar - Φbar_fortran), 
    abs(result1.Phibar - Φbar_fortran) / Φbar_fortran * 100)
println()

# 测试 2：使用 Julia 默认初始值
println("="^80)
println("测试 2：使用 Julia 默认初始值")
println("="^80)

result2 = solve(FixedMu(), T, μ)

println("\nJulia 求解结果（默认初始值）:")
@printf("  收敛状态: %s\n", result2.converged ? "✅ 收敛" : "❌ 未收敛")
@printf("  m_u     = %.6f fm⁻¹ = %.2f MeV\n", result2.masses[1], result2.masses[1] * ħc_MeV_fm)
@printf("  m_d     = %.6f fm⁻¹ = %.2f MeV\n", result2.masses[2], result2.masses[2] * ħc_MeV_fm)
@printf("  m_s     = %.6f fm⁻¹ = %.2f MeV\n", result2.masses[3], result2.masses[3] * ħc_MeV_fm)
@printf("  Φ       = %.8f\n", result2.Phi)
@printf("  Φ̄       = %.8f\n", result2.Phibar)
println()

# 对比两个 Julia 解
println("两个 Julia 解的差异:")
@printf("  Δm_u    = %.2e fm⁻¹ (%.4f%%)\n", 
    abs(result1.masses[1] - result2.masses[1]), 
    abs(result1.masses[1] - result2.masses[1]) / result2.masses[1] * 100)
@printf("  Δm_d    = %.2e fm⁻¹ (%.4f%%)\n", 
    abs(result1.masses[2] - result2.masses[2]), 
    abs(result1.masses[2] - result2.masses[2]) / result2.masses[2] * 100)
@printf("  Δm_s    = %.2e fm⁻¹ (%.4f%%)\n", 
    abs(result1.masses[3] - result2.masses[3]), 
    abs(result1.masses[3] - result2.masses[3]) / result2.masses[3] * 100)
@printf("  ΔΦ      = %.2e (%.4f%%)\n", 
    abs(result1.Phi - result2.Phi), 
    abs(result1.Phi - result2.Phi) / result2.Phi * 100)
@printf("  ΔΦ̄      = %.2e (%.4f%%)\n", 
    abs(result1.Phibar - result2.Phibar), 
    abs(result1.Phibar - result2.Phibar) / result2.Phibar * 100)
println()

# 测试 3：检查能隙方程的残差
println("="^80)
println("测试 3：检查能隙方程残差")
println("="^80)

# 使用 Fortran 的解计算残差
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))

using .GaussLegendre
using .OneLoopIntegrals

# 生成积分节点
nodes_p, weights_p = gauleg(0.0, 20.0, 32)

# 计算 A 函数
A_u = OneLoopIntegrals.A(m_u_fortran, μ, T, Φ_fortran, Φbar_fortran, nodes_p, weights_p)
A_d = OneLoopIntegrals.A(m_d_fortran, μ, T, Φ_fortran, Φbar_fortran, nodes_p, weights_p)
A_s = OneLoopIntegrals.A(m_s_fortran, μ, T, Φ_fortran, Φbar_fortran, nodes_p, weights_p)

println("\n使用 Fortran 解计算的 A 函数:")
@printf("  A_u     = %.6f fm⁻²\n", A_u)
@printf("  A_d     = %.6f fm⁻²\n", A_d)
@printf("  A_s     = %.6f fm⁻²\n", A_s)
println()

# 计算能隙方程残差
residual_u = m_u_fortran - m_ud0_inv_fm + 2.0 * G_fm2 * m_u_fortran * A_u - 
             4.0 * K_fm5 * m_u_fortran * A_u * A_d * A_s
residual_d = m_d_fortran - m_ud0_inv_fm + 2.0 * G_fm2 * m_d_fortran * A_d - 
             4.0 * K_fm5 * m_d_fortran * A_u * A_d * A_s
residual_s = m_s_fortran - m_s0_inv_fm + 2.0 * G_fm2 * m_s_fortran * A_s - 
             4.0 * K_fm5 * m_s_fortran * A_u * A_d * A_s

println("能隙方程残差（使用 Fortran 解）:")
@printf("  F_u     = %.2e\n", residual_u)
@printf("  F_d     = %.2e\n", residual_d)
@printf("  F_s     = %.2e\n", residual_s)
println()

# 对比 Fortran 的 A 值
println("与 Fortran A 值的对比:")
A_u_fortran = -4.972144
A_s_fortran = -5.165541

@printf("  A_u (Julia)   = %.6f fm⁻²\n", A_u)
@printf("  A_u (Fortran) = %.6f fm⁻²\n", A_u_fortran)
@printf("  差异          = %.2f%%\n", abs(A_u - A_u_fortran) / abs(A_u_fortran) * 100)
println()
@printf("  A_s (Julia)   = %.6f fm⁻²\n", A_s)
@printf("  A_s (Fortran) = %.6f fm⁻²\n", A_s_fortran)
@printf("  差异          = %.2f%%\n", abs(A_s - A_s_fortran) / abs(A_s_fortran) * 100)
println()

# 总结
println("="^80)
println("总结")
println("="^80)
println()

if abs(result1.masses[3] - m_s_fortran) / m_s_fortran < 0.01
    println("✅ 使用 Fortran 初始值，Julia 收敛到相同的解（差异 < 1%）")
    println("   → 说明 Julia 的求解器是正确的")
    println("   → 问题在于默认初始值导致收敛到不同的局部解")
elseif abs(result1.masses[3] - result2.masses[3]) / result2.masses[3] < 0.01
    println("❌ 无论使用什么初始值，Julia 都收敛到相同的解")
    println("   → 说明 Julia 的求解器总是收敛到同一个解")
    println("   → 问题可能在于能隙方程的实现不同")
else
    println("⚠️  Julia 收敛到不同的解，取决于初始值")
    println("   → 说明能隙方程可能有多个解")
    println("   → 需要检查哪个解是物理的")
end
println()

if abs(A_u - A_u_fortran) / abs(A_u_fortran) > 0.05
    println("❌ 即使使用相同的 m, Φ, Julia 的 A 函数也与 Fortran 不同")
    println("   → 说明 A 函数的实现可能不同")
    println("   → 需要检查积分方法、分布函数等")
else
    println("✅ 使用相同的 m, Φ, Julia 的 A 函数与 Fortran 一致")
    println("   → 说明 A 函数的实现是正确的")
end
println()

println("="^80)
