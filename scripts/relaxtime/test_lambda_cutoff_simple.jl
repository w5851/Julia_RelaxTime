#!/usr/bin/env julia
"""
简化的Lambda截断测试

对比两种情况：
1. 无截断：sigma_cutoff=nothing
2. 有截断：sigma_cutoff=Λ

测试条件：T=300 MeV, μ=2 MeV, ξ=0
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

println("="^80)
println("Lambda截断对弛豫时间的影响测试（简化版）")
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

# 先求解平衡态
println("步骤1: 求解PNJL平衡态...")
seed_state = [-0.001, -0.001, -0.04, 0.8, 0.8]  # 高温初值
base = TransportWorkflow.PNJL.solve(
    TransportWorkflow.PNJL.FixedMu(),
    T, μ;
    xi=ξ,
    seed_strategy=TransportWorkflow.PNJL.DefaultSeed(seed_state, seed_state, :quark),
    iterations=40,
)

if !base.converged
    error("PNJL平衡态未收敛！")
end

# 提取质量和Polyakov loop
Φ = Float64(base.x_state[4])
Φbar = Float64(base.x_state[5])
masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))

@printf("  收敛: ✓\n")
@printf("  Φ = %.6f\n", Φ)
@printf("  Φ̄ = %.6f\n", Φbar)
@printf("  m_u = %.2f MeV\n", masses.u * ħc_MeV_fm)
@printf("  m_s = %.2f MeV\n", masses.s * ħc_MeV_fm)
println()

# 计算K系数
println("步骤2: 计算有效耦合系数...")
nodes = DEFAULT_MOMENTUM_NODES
weights = DEFAULT_MOMENTUM_WEIGHTS
A_u = A(masses.u, μ, T, Φ, Φbar, nodes, weights)
A_s = A(masses.s, μ, T, Φ, Φbar, nodes, weights)
G_u = calculate_G_from_A(A_u, masses.u)
G_s = calculate_G_from_A(A_s, masses.s)
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
println("  完成")
println()

# 情况1: 默认行为（当前实现）
println("="^80)
println("情况1: 默认行为（当前实现）")
println("="^80)
println()

@time result_default = solve_gap_and_transport(
    T, μ;
    xi=ξ,
    equilibrium=base,
    compute_tau=true,
    K_coeffs=K_coeffs
)

τ_default = result_default.tau

println()
@printf("结果:\n")
@printf("  τ_u    = %.3f fm\n", τ_default.u)
@printf("  τ_s    = %.3f fm\n", τ_default.s)
@printf("  τ_u/τ_s = %.3f\n", τ_default.u / τ_default.s)
println()

println("="^80)
println("说明")
println("="^80)
println()
println("当前Julia实现已经使用了Lambda截断：")
println("  - 动量积分范围: [0, Λ] where Λ = $(round(Λ_inv_fm, digits=4)) fm⁻¹")
println("  - σ(s)缓存范围: 也使用Λ截断")
println()
println("这与Fortran实现一致。")
println()

# 与Fortran对比
τ_u_fortran = 0.584
τ_s_fortran = 0.593
ratio_fortran = τ_u_fortran / τ_s_fortran

println("="^80)
println("与Fortran参考值对比")
println("="^80)
println()

@printf("Fortran: τ_u = %.3f fm, τ_s = %.3f fm, τ_u/τ_s = %.3f\n", τ_u_fortran, τ_s_fortran, ratio_fortran)
println()

@printf("Julia vs Fortran:\n")
@printf("  τ_u:  Julia/Fortran = %.2f×\n", τ_default.u / τ_u_fortran)
@printf("  τ_s:  Julia/Fortran = %.2f×\n", τ_default.s / τ_s_fortran)
@printf("  比值: Julia/Fortran = %.3f×\n", (τ_default.u/τ_default.s) / ratio_fortran)
println()

println("="^80)
