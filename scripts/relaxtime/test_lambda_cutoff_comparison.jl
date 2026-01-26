#!/usr/bin/env julia
"""
测试Lambda截断对弛豫时间的影响

对比三种情况：
1. 原始实现：动量积分 [0, ∞)，无s截断
2. 新实现：动量积分 [0, ∞)，但σ(s)在s超出Lambda时归零
3. Fortran参考：动量积分 [0, Λ]

测试条件：T=300 MeV, μ=2 MeV, ξ=0
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))

using .Constants_PNJL: Λ_inv_fm, ħc_MeV_fm
using .TransportWorkflow: solve_gap_and_transport, FixedMu

println("="^80)
println("Lambda截断对弛豫时间的影响测试")
println("="^80)
println()

# 测试参数
T_MeV = 300.0
μ_MeV = 2.0
ξ = 0.0

T = T_MeV / ħc_MeV_fm
μ = μ_MeV / ħc_MeV_fm

println("测试条件:")
println("  T = $T_MeV MeV = $(round(T, digits=6)) fm⁻¹")
println("  μ = $μ_MeV MeV = $(round(μ, digits=6)) fm⁻¹")
println("  ξ = $ξ")
println("  Λ = $(round(Λ_inv_fm * ħc_MeV_fm, digits=2)) MeV = $(round(Λ_inv_fm, digits=4)) fm⁻¹")
println()

# 情况1: 原始实现（无截断）
println("="^80)
println("情况1: 动量积分 [0, ∞)，无s截断")
println("="^80)
println()

println("求解PNJL平衡态并计算弛豫时间...")
@time result_no_cut = solve_gap_and_transport(
    FixedMu(),
    T, μ;
    xi=ξ,
    compute_tau=true,
    tau_kwargs=(sigma_cutoff=nothing,)  # 不使用截断，通过tau_kwargs传递
)

τ_no_cut = result_no_cut.relaxation_times

println()
println("结果:")
println("  τ_u    = $(round(τ_no_cut.u, digits=3)) fm")
println("  τ_d    = $(round(τ_no_cut.d, digits=3)) fm")
println("  τ_s    = $(round(τ_no_cut.s, digits=3)) fm")
println("  τ_ubar = $(round(τ_no_cut.ubar, digits=3)) fm")
println("  τ_dbar = $(round(τ_no_cut.dbar, digits=3)) fm")
println("  τ_sbar = $(round(τ_no_cut.sbar, digits=3)) fm")
println()
println("  τ_u/τ_s = $(round(τ_no_cut.u / τ_no_cut.s, digits=3))")
println()

# 情况2: 新实现（σ(s)在s超出Lambda时归零）
println("="^80)
println("情况2: 动量积分 [0, ∞)，σ(s)在s超出Λ时归零")
println("="^80)
println()

println("求解PNJL平衡态并计算弛豫时间...")
@time result_with_cut = solve_gap_and_transport(
    FixedMu(),
    T, μ;
    xi=ξ,
    compute_tau=true,
    tau_kwargs=(sigma_cutoff=Λ_inv_fm,)  # 使用Lambda截断，通过tau_kwargs传递
)

τ_with_cut = result_with_cut.relaxation_times

println()
println("结果:")
println("  τ_u    = $(round(τ_with_cut.u, digits=3)) fm")
println("  τ_d    = $(round(τ_with_cut.d, digits=3)) fm")
println("  τ_s    = $(round(τ_with_cut.s, digits=3)) fm")
println("  τ_ubar = $(round(τ_with_cut.ubar, digits=3)) fm")
println("  τ_dbar = $(round(τ_with_cut.dbar, digits=3)) fm")
println("  τ_sbar = $(round(τ_with_cut.sbar, digits=3)) fm")
println()
println("  τ_u/τ_s = $(round(τ_with_cut.u / τ_with_cut.s, digits=3))")
println()

# 对比分析
println("="^80)
println("对比分析")
println("="^80)
println()

ratio_u = τ_with_cut.u / τ_no_cut.u
ratio_s = τ_with_cut.s / τ_no_cut.s

println("绝对值变化:")
println("  τ_u:  $(round(τ_no_cut.u, digits=3)) → $(round(τ_with_cut.u, digits=3)) fm (×$(round(ratio_u, digits=3)))")
println("  τ_s:  $(round(τ_no_cut.s, digits=3)) → $(round(τ_with_cut.s, digits=3)) fm (×$(round(ratio_s, digits=3)))")
println()

ratio_change = (τ_with_cut.u / τ_with_cut.s) / (τ_no_cut.u / τ_no_cut.s)
println("比值变化:")
println("  τ_u/τ_s: $(round(τ_no_cut.u / τ_no_cut.s, digits=3)) → $(round(τ_with_cut.u / τ_with_cut.s, digits=3)) (×$(round(ratio_change, digits=3)))")
println()

# 与Fortran对比
println("="^80)
println("与Fortran参考值对比")
println("="^80)
println()

τ_u_fortran = 0.584
τ_s_fortran = 0.593
ratio_fortran = τ_u_fortran / τ_s_fortran

println("Fortran参考值 (T=300 MeV, μ=2 MeV):")
println("  τ_u = $τ_u_fortran fm")
println("  τ_s = $τ_s_fortran fm")
println("  τ_u/τ_s = $(round(ratio_fortran, digits=3))")
println()

println("Julia (无截断) vs Fortran:")
println("  τ_u:  Julia/Fortran = $(round(τ_no_cut.u / τ_u_fortran, digits=2))×")
println("  τ_s:  Julia/Fortran = $(round(τ_no_cut.s / τ_s_fortran, digits=2))×")
println("  比值: Julia/Fortran = $(round((τ_no_cut.u/τ_no_cut.s) / ratio_fortran, digits=3))×")
println()

println("Julia (有截断) vs Fortran:")
println("  τ_u:  Julia/Fortran = $(round(τ_with_cut.u / τ_u_fortran, digits=2))×")
println("  τ_s:  Julia/Fortran = $(round(τ_with_cut.s / τ_s_fortran, digits=2))×")
println("  比值: Julia/Fortran = $(round((τ_with_cut.u/τ_with_cut.s) / ratio_fortran, digits=3))×")
println()

println("="^80)
println("结论")
println("="^80)
println()

if abs(ratio_u - 1.0) > 0.1 || abs(ratio_s - 1.0) > 0.1
    println("✓ Lambda截断对弛豫时间有显著影响")
    println("  - 绝对值变化: $(round((1-ratio_u)*100, digits=1))% (u), $(round((1-ratio_s)*100, digits=1))% (s)")
else
    println("⚠ Lambda截断对弛豫时间影响较小")
end

if abs(ratio_change - 1.0) > 0.05
    println("✓ Lambda截断改变了τ_u/τ_s比值")
    println("  - 比值变化: $(round((ratio_change-1)*100, digits=1))%")
else
    println("⚠ Lambda截断对τ_u/τ_s比值影响较小")
end

println()
println("与Fortran的差异:")
julia_fortran_ratio_no_cut = (τ_no_cut.u/τ_no_cut.s) / ratio_fortran
julia_fortran_ratio_with_cut = (τ_with_cut.u/τ_with_cut.s) / ratio_fortran

if abs(julia_fortran_ratio_with_cut - 1.0) < abs(julia_fortran_ratio_no_cut - 1.0)
    println("✓ 使用Lambda截断后，Julia与Fortran的比值更接近")
    println("  - 改进: $(round(abs(julia_fortran_ratio_no_cut - 1.0)*100, digits=1))% → $(round(abs(julia_fortran_ratio_with_cut - 1.0)*100, digits=1))%")
else
    println("⚠ Lambda截断未能改善与Fortran的一致性")
end

println()
println("="^80)
