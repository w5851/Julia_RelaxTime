#!/usr/bin/env julia
"""
调试 ssbar_to_uubar 过程

这个过程的散射率比Fortran大3.238倍，需要找出原因
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "AverageScatteringRate.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))

using .Constants_PNJL: Λ_inv_fm, ħc_MeV_fm, G_fm2, K_fm5, SCATTERING_MESON_MAP
using .TransportWorkflow: solve_gap_and_transport
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .EffectiveCouplings.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .AverageScatteringRate: average_scattering_rate, CrossSectionCache, build_w0cdf_pchip_cache
using .TotalCrossSection: total_cross_section
using Printf

println("="^80)
println("调试 ssbar_to_uubar 过程")
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

# 检查过程定义
println("="^80)
println("过程定义")
println("="^80)
println()

process = :ssbar_to_uubar
process_info = SCATTERING_MESON_MAP[process]

println("过程: ssbar_to_uubar (s + s̄ → u + ū)")
println("类型: ", process_info[:type])
println()

println("通道:")
for (channel, info) in process_info[:channels]
    println("  $channel 通道:")
    println("    simple介子: ", info[:simple])
    println("    mixed_P: ", info[:mixed_P])
    println("    mixed_S: ", info[:mixed_S])
end
println()

# 对比过程：uubar_to_ssbar
println("对比过程: uubar_to_ssbar (u + ū → s + s̄)")
process_comp = :uubar_to_ssbar
process_comp_info = SCATTERING_MESON_MAP[process_comp]
println("类型: ", process_comp_info[:type])
println()

println("通道:")
for (channel, info) in process_comp_info[:channels]
    println("  $channel 通道:")
    println("    simple介子: ", info[:simple])
    println("    mixed_P: ", info[:mixed_P])
    println("    mixed_S: ", info[:mixed_S])
end
println()

# 计算散射率
println("="^80)
println("散射率计算")
println("="^80)
println()

println("计算 ssbar_to_uubar...")
w_ssbar_uubar = average_scattering_rate(
    :ssbar_to_uubar,
    quark_params,
    thermo_params,
    K_coeffs;
    p_nodes=20,
    angle_nodes=4,
    phi_nodes=8,
    sigma_cutoff=Λ_inv_fm
)

println("计算 uubar_to_ssbar...")
w_uubar_ssbar = average_scattering_rate(
    :uubar_to_ssbar,
    quark_params,
    thermo_params,
    K_coeffs;
    p_nodes=20,
    angle_nodes=4,
    phi_nodes=8,
    sigma_cutoff=Λ_inv_fm
)

println()
@printf("Julia结果:\n")
@printf("  w̄(ssbar→uubar) = %.6e fm⁻³\n", w_ssbar_uubar)
@printf("  w̄(uubar→ssbar) = %.6e fm⁻³\n", w_uubar_ssbar)
@printf("  比值 = %.3f\n", w_ssbar_uubar / w_uubar_ssbar)
println()

@printf("Fortran结果:\n")
@printf("  w̄(ssbar→uubar) = %.6e fm⁻³\n", 0.6074054988e-01)
@printf("  w̄(uubar→ssbar) = %.6e fm⁻³\n", 0.4692409224e-01)
@printf("  比值 = %.3f\n", 0.6074054988e-01 / 0.4692409224e-01)
println()

@printf("Julia/Fortran:\n")
@printf("  ssbar→uubar: %.3f\n", w_ssbar_uubar / 0.6074054988e-01)
@printf("  uubar→ssbar: %.3f\n", w_uubar_ssbar / 0.4692409224e-01)
println()

# 检查详细平衡
println("="^80)
println("详细平衡检查")
println("="^80)
println()

println("详细平衡原理:")
println("  w̄(i→f) / w̄(f→i) = exp(-(E_f - E_i)/T)")
println("  对于 ssbar ↔ uubar:")
println("    E_ssbar ≈ 2*m_s = ", 2*masses.s*ħc_MeV_fm, " MeV")
println("    E_uubar ≈ 2*m_u = ", 2*masses.u*ħc_MeV_fm, " MeV")
println("    ΔE = E_uubar - E_ssbar = ", (2*masses.u - 2*masses.s)*ħc_MeV_fm, " MeV")
println("    exp(-ΔE/T) = ", exp(-(2*masses.u - 2*masses.s)/T))
println()

@printf("Julia比值:\n")
@printf("  w̄(uubar→ssbar) / w̄(ssbar→uubar) = %.6f\n", w_uubar_ssbar / w_ssbar_uubar)
@printf("  理论值 exp(-ΔE/T) = %.6f\n", exp(-(2*masses.u - 2*masses.s)/T))
println()

@printf("Fortran比值:\n")
@printf("  w̄(uubar→ssbar) / w̄(ssbar→uubar) = %.6f\n", 0.4692409224e-01 / 0.6074054988e-01)
@printf("  理论值 exp(-ΔE/T) = %.6f\n", exp(-(2*masses.u - 2*masses.s)/T))
println()

# 计算截面
println("="^80)
println("截面对比")
println("="^80)
println()

# 在阈值附近计算截面
s_threshold_ssbar = (masses.s + masses.s)^2
s_threshold_uubar = (masses.u + masses.u)^2

println("阈值:")
@printf("  s_threshold(ssbar) = %.6f fm⁻² = %.2f MeV²\n", s_threshold_ssbar, s_threshold_ssbar * ħc_MeV_fm^2)
@printf("  s_threshold(uubar) = %.6f fm⁻² = %.2f MeV²\n", s_threshold_uubar, s_threshold_uubar * ħc_MeV_fm^2)
println()

# 在几个s值处计算截面
s_values = [s_threshold_ssbar * 1.1, s_threshold_ssbar * 2.0, s_threshold_ssbar * 5.0]

println("截面 σ(s):")
println("-"^80)
@printf("%-20s %15s %15s %12s\n", "s (MeV²)", "σ(ssbar→uubar)", "σ(uubar→ssbar)", "比值")
println("-"^80)

for s in s_values
    σ_ssbar_uubar = total_cross_section(:ssbar_to_uubar, s, quark_params, thermo_params, K_coeffs)
    σ_uubar_ssbar = total_cross_section(:uubar_to_ssbar, s, quark_params, thermo_params, K_coeffs)
    ratio = σ_ssbar_uubar / σ_uubar_ssbar
    @printf("%20.2f %15.6e %15.6e %12.3f\n", s * ħc_MeV_fm^2, σ_ssbar_uubar, σ_uubar_ssbar, ratio)
end

println()
println("="^80)
println("结论")
println("="^80)
println()

if w_ssbar_uubar / 0.6074054988e-01 > 2.0
    println("⚠️ Julia的 ssbar→uubar 散射率比Fortran大 $(round(w_ssbar_uubar / 0.6074054988e-01, digits=2))倍")
    println()
    println("可能的原因:")
    println("  1. 过程定义错误（粒子映射）")
    println("  2. 介子通道定义错误")
    println("  3. 散射振幅计算错误")
    println("  4. 积分方法差异")
    println()
else
    println("✓ Julia和Fortran的散射率基本一致")
    println()
end

println("="^80)
