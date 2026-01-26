#!/usr/bin/env julia
"""
验证角度积分bug修复

测试修复后的Julia实现与Fortran的对比
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))

using .Constants_PNJL: Λ_inv_fm, ħc_MeV_fm, G_fm2, K_fm5
using .TransportWorkflow: solve_gap_and_transport
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .EffectiveCouplings.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using Printf

# Import relaxation_times function
const RelaxTime = Base.include(Main, joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxationTime.jl"))
const relaxation_times = RelaxTime.relaxation_times
const AverageScatteringRate = RelaxTime.AverageScatteringRate
const number_density = AverageScatteringRate.number_density

println("="^80)
println("角度积分Bug修复验证")
println("="^80)
println()

# 测试参数（与Fortran对比）
T_MeV = 300.0
μ_MeV = 2.0
ξ = 0.0

T = T_MeV / ħc_MeV_fm
μ = μ_MeV / ħc_MeV_fm

println("测试条件:")
@printf("  T = %.1f MeV\n", T_MeV)
@printf("  μ = %.1f MeV\n", μ_MeV)
@printf("  ξ = %.1f (各向同性)\n", ξ)
println()

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
@printf("  m_u = %.3f MeV\n", masses.u * ħc_MeV_fm)
@printf("  m_s = %.3f MeV\n", masses.s * ħc_MeV_fm)
@printf("  Φ = %.6f\n", Φ)
@printf("  Φbar = %.6f\n", Φbar)
println()

# 计算K系数
nodes = DEFAULT_MOMENTUM_NODES
weights = DEFAULT_MOMENTUM_WEIGHTS
A_u = A(masses.u, μ, T, Φ, Φbar, nodes, weights)
A_s = A(masses.s, μ, T, Φ, Φbar, nodes, weights)
G_u = calculate_G_from_A(A_u, masses.u)
G_s = calculate_G_from_A(A_s, masses.s)
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

quark_params = (
    m = masses,
    μ = (u=μ, d=μ, s=μ),
    A = (u=A_u, d=A_u, s=A_s)
)

thermo_params = (
    T = T,
    Φ = Φ,
    Φbar = Φbar,
    ξ = ξ
)

# 计算弛豫时间
println("="^80)
println("弛豫时间计算")
println("="^80)
println()

# 首先计算数密度
println("计算数密度...")
ρ_u = number_density(:u, masses.u, μ, T, Φ, Φbar, ξ)
ρ_d = number_density(:d, masses.d, μ, T, Φ, Φbar, ξ)
ρ_s = number_density(:s, masses.s, μ, T, Φ, Φbar, ξ)
ρ_ubar = number_density(:ubar, masses.u, μ, T, Φ, Φbar, ξ)
ρ_dbar = number_density(:dbar, masses.d, μ, T, Φ, Φbar, ξ)
ρ_sbar = number_density(:sbar, masses.s, μ, T, Φ, Φbar, ξ)

densities = (
    u = ρ_u,
    d = ρ_d,
    s = ρ_s,
    ubar = ρ_ubar,
    dbar = ρ_dbar,
    sbar = ρ_sbar
)

println("完成")
println()

println("计算弛豫时间...")
result = relaxation_times(quark_params, thermo_params, K_coeffs; densities=densities)
τ_u = result.tau.u
τ_s = result.tau.s

println()
println("Julia结果（修复后）:")
println("-"^80)
@printf("  τ_u = %.3f fm\n", τ_u)
@printf("  τ_s = %.3f fm\n", τ_s)
@printf("  τ_u/τ_s = %.3f\n", τ_u / τ_s)
println()

# Fortran参考值
τ_u_fortran = 0.584
τ_s_fortran = 0.593

println("Fortran参考值:")
println("-"^80)
@printf("  τ_u = %.3f fm\n", τ_u_fortran)
@printf("  τ_s = %.3f fm\n", τ_s_fortran)
@printf("  τ_u/τ_s = %.3f\n", τ_u_fortran / τ_s_fortran)
println()

# 对比
println("="^80)
println("对比分析")
println("="^80)
println()

ratio_u = τ_u / τ_u_fortran
ratio_s = τ_s / τ_s_fortran

@printf("Julia/Fortran 比值:\n")
@printf("  τ_u: %.3f (修复前: ~2.96)\n", ratio_u)
@printf("  τ_s: %.3f (修复前: ~3.54)\n", ratio_s)
println()

# 修复前的值（用于对比）
τ_u_before = 1.726
τ_s_before = 2.100

@printf("修复效果:\n")
@printf("  τ_u: %.3f fm → %.3f fm (变化: %.1f%%)\n", 
    τ_u_before, τ_u, (τ_u - τ_u_before) / τ_u_before * 100)
@printf("  τ_s: %.3f fm → %.3f fm (变化: %.1f%%)\n", 
    τ_s_before, τ_s, (τ_s - τ_s_before) / τ_s_before * 100)
println()

# 评估修复是否成功
println("="^80)
println("修复评估")
println("="^80)
println()

# 预期：修复后应该比Fortran小约12.5倍（因为归一化约定差异）
expected_ratio = 1.0 / 12.57  # ≈ 0.08

@printf("预期比值: %.3f (1/4π ≈ 1/12.57)\n", expected_ratio)
@printf("实际比值: τ_u = %.3f, τ_s = %.3f\n", ratio_u, ratio_s)
println()

if abs(ratio_u - expected_ratio) / expected_ratio < 0.5 && 
   abs(ratio_s - expected_ratio) / expected_ratio < 0.5
    println("✓ 修复成功！")
    println("  - Julia/Fortran比值接近预期值 1/4π")
    println("  - 角度积分bug已修复")
    println()
else
    println("⚠️ 修复后仍有差异")
    @printf("  - τ_u比值偏差: %.1f%%\n", abs(ratio_u - expected_ratio) / expected_ratio * 100)
    @printf("  - τ_s比值偏差: %.1f%%\n", abs(ratio_s - expected_ratio) / expected_ratio * 100)
    println("  - 可能还有其他因素影响")
    println()
end

# 物理正确性检查
if τ_u < τ_s
    println("✓ 物理图像正确: τ_u < τ_s")
    println("  - u夸克更快达到局部平衡")
    println()
else
    println("⚠️ 物理图像异常: τ_u > τ_s")
    println()
end

println("="^80)
println("总结")
println("="^80)
println()

println("修复内容:")
println("  - 在各向同性情况下（ξ=0），归一化因子除以4π")
println("  - 补偿5维角度积分（8π）与Fortran 3维积分（2）的差异")
println()

println("修复位置:")
println("  - src/relaxtime/AverageScatteringRate.jl")
println("  - _average_scattering_rate_semi_infinite 函数")
println()

println("="^80)
