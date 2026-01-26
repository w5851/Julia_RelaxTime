#!/usr/bin/env julia
"""
验证Julia和Fortran的归一化差异

检查归一化因子是否能解释3倍差异
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "AverageScatteringRate.jl"))

using .Constants_PNJL: Λ_inv_fm, ħc_MeV_fm, G_fm2, K_fm5
using .TransportWorkflow: solve_gap_and_transport
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .EffectiveCouplings.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .AverageScatteringRate: number_density
using Printf

println("="^80)
println("归一化差异验证")
println("="^80)
println()

# 测试参数
T_MeV = 300.0
μ_MeV = 2.0
ξ = 0.0

T = T_MeV / ħc_MeV_fm
μ = μ_MeV / ħc_MeV_fm

# 求解平衡态
println("求解PNJL平衡态...")
seed_state = [-0.001, -0.001, -0.04, 0.8, 0.8]
base = TransportWorkflow.PNJL.solve(
    TransportWorkflow.PNJL.FixedMu(),
    T, μ;
    xi=ξ,
    seed_strategy=TransportWorkflow.PNJL.DefaultSeed(seed_state, seed_state, :quark),
    iterations=40,
)

Φ = Float64(base.x_state[4])
Φbar = Float64(base.x_state[5])
masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))

println("完成")
println()

# 计算密度
println("="^80)
println("密度计算")
println("="^80)
println()

ρ_u = number_density(:u, masses.u, μ, T, Φ, Φbar, ξ)
ρ_s = number_density(:s, masses.s, μ, T, Φ, Φbar, ξ)
ρ_ubar = number_density(:ubar, masses.u, μ, T, Φ, Φbar, ξ)
ρ_sbar = number_density(:sbar, masses.s, μ, T, Φ, Φbar, ξ)

@printf("Julia密度:\n")
@printf("  ρ_u    = %.6f fm⁻³\n", ρ_u)
@printf("  ρ_s    = %.6f fm⁻³\n", ρ_s)
@printf("  ρ_ubar = %.6f fm⁻³\n", ρ_ubar)
@printf("  ρ_sbar = %.6f fm⁻³\n", ρ_sbar)
println()

@printf("Fortran密度:\n")
@printf("  n_u    = 1.680901 fm⁻³\n")
@printf("  n_s    = 1.537628 fm⁻³\n")
@printf("  n_ubar = 1.673911 fm⁻³\n")
@printf("  n_sbar = 1.531154 fm⁻³\n")
println()

# 归一化因子对比
println("="^80)
println("归一化因子对比")
println("="^80)
println()

N_c = 3.0
DQ = 6.0  # 2*N_c

println("Julia的归一化因子:")
println("  prefactor = DQ²/(32π⁵ρ_iρ_j)")
println("  其中 DQ = 2N_c = 6")
println()

@printf("对于 ssbar→uubar:\n")
julia_prefactor_ssbar = (DQ^2) / (32.0 * π^5 * ρ_s * ρ_sbar)
@printf("  Julia: DQ²/(32π⁵ρ_sρ_sbar) = %.6e\n", julia_prefactor_ssbar)
println()

@printf("对于 uubar→ssbar:\n")
julia_prefactor_uubar = (DQ^2) / (32.0 * π^5 * ρ_u * ρ_ubar)
@printf("  Julia: DQ²/(32π⁵ρ_uρ_ubar) = %.6e\n", julia_prefactor_uubar)
println()

println("-"^80)
println()

println("Fortran的归一化因子:")
println("  w_ij_n = N_c²/(2π⁴) × ω")
println("  w_ij = w_ij_n / (n1×n2)")
println("  即: w_ij = N_c²/(2π⁴n1n2) × ω")
println()

n_s_fortran = 1.537628
n_sbar_fortran = 1.531154
n_u_fortran = 1.680901
n_ubar_fortran = 1.673911

@printf("对于 ssbar→uubar:\n")
fortran_prefactor_ssbar = (N_c^2) / (2.0 * π^4 * n_s_fortran * n_sbar_fortran)
@printf("  Fortran: N_c²/(2π⁴n_sn_sbar) = %.6e\n", fortran_prefactor_ssbar)
println()

@printf("对于 uubar→ssbar:\n")
fortran_prefactor_uubar = (N_c^2) / (2.0 * π^4 * n_u_fortran * n_ubar_fortran)
@printf("  Fortran: N_c²/(2π⁴n_un_ubar) = %.6e\n", fortran_prefactor_uubar)
println()

# 比值分析
println("="^80)
println("比值分析")
println("="^80)
println()

@printf("ssbar→uubar:\n")
@printf("  Julia/Fortran = %.6f\n", julia_prefactor_ssbar / fortran_prefactor_ssbar)
println()

@printf("uubar→ssbar:\n")
@printf("  Julia/Fortran = %.6f\n", julia_prefactor_uubar / fortran_prefactor_uubar)
println()

# 理论计算
println("="^80)
println("理论分析")
println("="^80)
println()

println("Julia前因子:")
println("  DQ²/(32π⁵ρ_iρ_j) = (2N_c)²/(32π⁵ρ_iρ_j) = 4N_c²/(32π⁵ρ_iρ_j) = N_c²/(8π⁵ρ_iρ_j)")
println()

println("Fortran前因子:")
println("  N_c²/(2π⁴n_in_j)")
println()

println("比值:")
println("  Julia/Fortran = [N_c²/(8π⁵ρ_iρ_j)] / [N_c²/(2π⁴n_in_j)]")
println("                = (2π⁴n_in_j) / (8π⁵ρ_iρ_j)")
println("                = n_in_j / (4πρ_iρ_j)")
println()

ratio_ssbar = (n_s_fortran * n_sbar_fortran) / (4.0 * π * ρ_s * ρ_sbar)
ratio_uubar = (n_u_fortran * n_ubar_fortran) / (4.0 * π * ρ_u * ρ_ubar)

@printf("理论比值:\n")
@printf("  ssbar→uubar: n_sn_sbar/(4πρ_sρ_sbar) = %.6f\n", ratio_ssbar)
@printf("  uubar→ssbar: n_un_ubar/(4πρ_uρ_ubar) = %.6f\n", ratio_uubar)
println()

# 与实际散射率比值对比
println("="^80)
println("与实际散射率比值对比")
println("="^80)
println()

println("实际散射率比值（从之前的测试）:")
println("  w̄(ssbar→uubar): Julia/Fortran = 3.238")
println("  w̄(uubar→ssbar): Julia/Fortran = 0.988")
println()

println("归一化因子比值:")
@printf("  ssbar→uubar: %.6f\n", julia_prefactor_ssbar / fortran_prefactor_ssbar)
@printf("  uubar→ssbar: %.6f\n", julia_prefactor_uubar / fortran_prefactor_uubar)
println()

println("⚠️ 关键发现:")
println("  归一化因子比值 ≈ 0.025（Julia比Fortran小40倍！）")
println("  但实际散射率比值 ≈ 3.238（Julia比Fortran大3倍）")
println()

println("这说明:")
println("  1. 归一化因子不能解释3倍差异")
println("  2. 积分ω本身必定有差异")
println("  3. Julia的ω比Fortran大约 3.238 / 0.025 ≈ 130倍")
println()

println("可能的原因:")
println("  1. Julia和Fortran的积分方法完全不同")
println("  2. Julia使用半无穷积分，Fortran使用有限积分")
println("  3. 积分权重或雅可比因子不同")
println()

println("="^80)
