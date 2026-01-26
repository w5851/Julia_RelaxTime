#!/usr/bin/env julia
"""
测试t积分的收敛性

检查 ssbar_to_uubar 和 uubar_to_ssbar 的t积分是否收敛
特别关注6个节点是否足够
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))

using .Constants_PNJL: Λ_inv_fm, ħc_MeV_fm, G_fm2, K_fm5
using .TransportWorkflow: solve_gap_and_transport
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .EffectiveCouplings.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .TotalCrossSection: total_cross_section
using Printf

println("="^80)
println("t积分收敛性测试")
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

println("完成")
println()

# 测试不同的t积分节点数
n_points_list = [4, 6, 8, 12, 16, 24, 32, 48, 64]

# 测试s值：在ssbar阈值附近
s_threshold_ssbar = (masses.s + masses.s)^2
s_test = s_threshold_ssbar * 1.1

println("="^80)
println("收敛性测试")
println("="^80)
println()

@printf("测试条件:\n")
@printf("  s = %.2f MeV² (%.2f × s_threshold)\n", s_test * ħc_MeV_fm^2, s_test / s_threshold_ssbar)
println()

# ssbar -> uubar
println("过程1: ssbar → uubar")
println("-"^80)
@printf("%-10s %20s %15s %15s\n", "n_points", "σ (fm²)", "相对变化", "计算时间")
println("-"^80)

σ_values_1 = Float64[]
for (i, n) in enumerate(n_points_list)
    t_start = time()
    σ = total_cross_section(:ssbar_to_uubar, s_test, quark_params, thermo_params, K_coeffs; n_points=n)
    t_elapsed = time() - t_start
    
    push!(σ_values_1, σ)
    rel_change = i > 1 ? abs(σ - σ_values_1[i-1]) / σ_values_1[i-1] * 100 : 0.0
    @printf("%-10d %20.6e %14.2f%% %14.3f ms\n", n, σ, rel_change, t_elapsed * 1000)
end
println()

# uubar -> ssbar
println("过程2: uubar → ssbar")
println("-"^80)
@printf("%-10s %20s %15s %15s\n", "n_points", "σ (fm²)", "相对变化", "计算时间")
println("-"^80)

σ_values_2 = Float64[]
for (i, n) in enumerate(n_points_list)
    t_start = time()
    σ = total_cross_section(:uubar_to_ssbar, s_test, quark_params, thermo_params, K_coeffs; n_points=n)
    t_elapsed = time() - t_start
    
    push!(σ_values_2, σ)
    rel_change = i > 1 ? abs(σ - σ_values_2[i-1]) / σ_values_2[i-1] * 100 : 0.0
    @printf("%-10d %20.6e %14.2f%% %14.3f ms\n", n, σ, rel_change, t_elapsed * 1000)
end
println()

# 对比比值
println("="^80)
println("比值分析")
println("="^80)
println()

println("σ(ssbar→uubar) / σ(uubar→ssbar) 随节点数的变化:")
println("-"^80)
@printf("%-10s %20s %20s %15s\n", "n_points", "σ(ssbar→uubar)", "σ(uubar→ssbar)", "比值")
println("-"^80)

for n in n_points_list
    σ_1 = total_cross_section(:ssbar_to_uubar, s_test, quark_params, thermo_params, K_coeffs; n_points=n)
    σ_2 = total_cross_section(:uubar_to_ssbar, s_test, quark_params, thermo_params, K_coeffs; n_points=n)
    ratio = σ_1 / σ_2
    @printf("%-10d %20.6e %20.6e %15.3f\n", n, σ_1, σ_2, ratio)
end
println()

# 检查Fortran使用的节点数
println("="^80)
println("与Fortran对比")
println("="^80)
println()

println("Fortran可能使用的t积分节点数：")
println("  - 需要检查Fortran代码中的积分节点数")
println("  - 如果Fortran使用更多节点，可能导致更准确的结果")
println()

# 推荐的节点数
println("="^80)
println("结论")
println("="^80)
println()

# 计算收敛性
σ_6 = total_cross_section(:ssbar_to_uubar, s_test, quark_params, thermo_params, K_coeffs; n_points=6)
σ_64 = total_cross_section(:ssbar_to_uubar, s_test, quark_params, thermo_params, K_coeffs; n_points=64)
convergence_error = abs(σ_6 - σ_64) / σ_64 * 100

@printf("收敛性分析 (ssbar→uubar):\n")
@printf("  σ(n=6)  = %.6e fm²\n", σ_6)
@printf("  σ(n=64) = %.6e fm²\n", σ_64)
@printf("  相对误差 = %.2f%%\n", convergence_error)
println()

if convergence_error > 5.0
    println("⚠️ 警告：6个节点的积分未收敛！")
    println("   相对误差 > 5%，建议增加节点数")
    println()
    
    # 找到收敛的节点数
    for n in [8, 12, 16, 24, 32]
        σ_n = total_cross_section(:ssbar_to_uubar, s_test, quark_params, thermo_params, K_coeffs; n_points=n)
        err = abs(σ_n - σ_64) / σ_64 * 100
        if err < 1.0
            @printf("   推荐使用 n = %d 个节点（误差 %.2f%%）\n", n, err)
            break
        end
    end
    println()
    
    println("这可能解释了与Fortran的3倍差异！")
    println("  - 如果Fortran使用更多节点，积分更准确")
    println("  - Julia使用6个节点可能不够")
    println()
else
    println("✓ 6个节点的积分已收敛（误差 < 5%）")
    println()
    println("积分收敛性不是问题的根源")
    println("  - 需要检查其他因素（归一化、权重等）")
    println()
end

println("="^80)
