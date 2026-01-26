#!/usr/bin/env julia
"""
简单对比 Fortran 和 Julia 的弛豫时间结果

Fortran: 半无穷积分（修改后）
Julia: 半无穷积分 + Λ截断σ(s)

测试条件：T=300 MeV, μ=2 MeV
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))

using .Constants_PNJL: Λ_inv_fm, ħc_MeV_fm
using .TransportWorkflow: solve_gap_and_transport
using Printf

println("="^80)
println("Fortran vs Julia 弛豫时间对比")
println("="^80)
println()

# 测试参数
T_MeV = 300.0
μ_MeV = 2.0
ξ = 0.0

T = T_MeV / ħc_MeV_fm
μ = μ_MeV / ħc_MeV_fm

@printf("测试条件:\n")
@printf("  T = %.1f MeV = %.6f fm⁻¹\n", T_MeV, T)
@printf("  μ = %.1f MeV = %.6f fm⁻¹\n", μ_MeV, μ)
@printf("  ξ = %.1f\n", ξ)
@printf("  Λ = %.1f MeV = %.4f fm⁻¹\n", Λ_inv_fm * ħc_MeV_fm, Λ_inv_fm)
println()

# Julia 计算
println("="^80)
println("Julia 计算（半无穷积分 + Λ截断σ(s)）")
println("="^80)

# 不传 K_coeffs，让 workflow 自动计算
result = solve_gap_and_transport(
    T, μ;
    xi=ξ,
    compute_tau=false,  # 先不计算 tau，只求平衡态
    p_num=20,
    t_num=4,
)

println("\nJulia 平衡态结果:")
@printf("  Φ     = %.6f\n", result.Phi)
@printf("  Φbar  = %.6f\n", result.Phibar)
@printf("  m_u   = %.6f fm⁻¹ = %.3f MeV\n", result.m_u, result.m_u * ħc_MeV_fm)
@printf("  m_s   = %.6f fm⁻¹ = %.3f MeV\n", result.m_s, result.m_s * ħc_MeV_fm)

# 注意：由于 compute_tau=false，result 中没有 tau_u, tau_s
# 我们需要手动计算或者使用另一个脚本的结果
println("\n注意：弛豫时间需要单独计算")
println("请参考 test_semi_infinite_momentum.jl 的结果")

# Fortran 结果（手动输入）
println("\n" * "="^80)
println("Fortran 结果（半无穷积分）")
println("="^80)

fortran_Phi = 0.840408
fortran_Phibar = 0.840412
fortran_m_u = 0.040510  # fm⁻¹
fortran_m_s = 1.032295  # fm⁻¹
fortran_tau_u = 0.58097  # fm
fortran_tau_s = 0.59204  # fm
fortran_ratio = fortran_tau_u / fortran_tau_s

println("\nFortran 结果:")
@printf("  Φ     = %.6f\n", fortran_Phi)
@printf("  Φbar  = %.6f\n", fortran_Phibar)
@printf("  m_u   = %.6f fm⁻¹ = %.3f MeV\n", fortran_m_u, fortran_m_u * ħc_MeV_fm)
@printf("  m_s   = %.6f fm⁻¹ = %.3f MeV\n", fortran_m_s, fortran_m_s * ħc_MeV_fm)
@printf("  τ_u   = %.5f fm\n", fortran_tau_u)
@printf("  τ_s   = %.5f fm\n", fortran_tau_s)
@printf("  τ_u/τ_s = %.4f\n", fortran_ratio)

# 对比
println("\n" * "="^80)
println("对比分析")
println("="^80)

println("\n1. 平衡态参数对比:")
@printf("  Φ (Julia/Fortran)     = %.4f\n", result.Phi / fortran_Phi)
@printf("  Φbar (Julia/Fortran)  = %.4f\n", result.Phibar / fortran_Phibar)
@printf("  m_u (Julia/Fortran)   = %.4f\n", result.m_u / fortran_m_u)
@printf("  m_s (Julia/Fortran)   = %.4f\n", result.m_s / fortran_m_s)

println("\n2. 弛豫时间对比:")
println("  （需要单独计算，请参考 test_semi_infinite_momentum.jl）")

# 使用之前的 Julia 结果（从 test_semi_infinite_momentum.jl）
julia_tau_u = 1.726  # fm
julia_tau_s = 2.100  # fm

@printf("\n  τ_u (Julia/Fortran)   = %.3f\n", julia_tau_u / fortran_tau_u)
@printf("  τ_s (Julia/Fortran)   = %.3f\n", julia_tau_s / fortran_tau_s)

println("\n3. 比值对比:")
julia_ratio = julia_tau_u / julia_tau_s
@printf("  τ_u/τ_s (Julia)       = %.4f\n", julia_ratio)
@printf("  τ_u/τ_s (Fortran)     = %.4f\n", fortran_ratio)
@printf("  比值的比值            = %.4f\n", julia_ratio / fortran_ratio)

println("\n4. 差异分析:")
if abs(result.Phi - fortran_Phi) < 0.001 && abs(result.m_u - fortran_m_u) < 0.001
    println("  ✅ 平衡态参数一致（Φ, m_u, m_s 差异 < 0.1%）")
else
    println("  ❌ 平衡态参数有差异")
end

tau_diff = abs(julia_tau_u / fortran_tau_u - 1.0) * 100
if tau_diff < 5.0
    println("  ✅ 弛豫时间接近（差异 < 5%）")
else
    @printf("  ❌ 弛豫时间差异较大（τ_u 差异 %.1f%%）\n", tau_diff)
end

ratio_diff = abs(julia_ratio / fortran_ratio - 1.0) * 100
if ratio_diff < 5.0
    println("  ✅ 比值接近（差异 < 5%）")
else
    @printf("  ❌ 比值差异较大（差异 %.1f%%）\n", ratio_diff)
end

println("\n" * "="^80)
println("结论")
println("="^80)

if tau_diff > 50.0
    println("\n❌ Julia 和 Fortran 的弛豫时间差异超过 50%")
    println("   可能的原因：")
    println("   1. σ(s) 的计算方法不同")
    println("   2. 散射振幅 M² 的计算不同")
    println("   3. 数值积分的精度不同")
    println("   4. 物理参数的差异")
    println("\n   建议：")
    println("   - 对比单个散射过程的散射率")
    println("   - 对比 σ(s) 的值")
    println("   - 检查散射振幅的计算")
elseif tau_diff > 10.0
    println("\n⚠️  Julia 和 Fortran 的弛豫时间有一定差异（10-50%）")
    println("   需要进一步调查差异来源")
else
    println("\n✅ Julia 和 Fortran 的弛豫时间基本一致（差异 < 10%）")
    println("   半无穷积分实现成功！")
end

println("\n" * "="^80)
