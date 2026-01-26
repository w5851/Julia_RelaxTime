#!/usr/bin/env julia
"""
对比 ssbar_to_uubar 和 uubar_to_ssbar 的散射振幅

目标：找出为什么 ssbar_to_uubar 比 uubar_to_ssbar 大3倍
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "ScatteringAmplitude.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))

using .Constants_PNJL: Λ_inv_fm, ħc_MeV_fm, G_fm2, K_fm5, SCATTERING_MESON_MAP
using .TransportWorkflow: solve_gap_and_transport
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .EffectiveCouplings.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .ScatteringAmplitude: scattering_amplitude_squared
using .TotalCrossSection: total_cross_section
using Printf

println("="^80)
println("ssbar_to_uubar vs uubar_to_ssbar 散射振幅对比")
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
println("过程定义对比")
println("="^80)
println()

proc1 = :ssbar_to_uubar
proc2 = :uubar_to_ssbar

info1 = SCATTERING_MESON_MAP[proc1]
info2 = SCATTERING_MESON_MAP[proc2]

println("ssbar_to_uubar:")
println("  类型: ", info1[:type])
for (ch, chinfo) in info1[:channels]
    println("  $ch 通道:")
    println("    simple: ", chinfo[:simple])
    println("    mixed_P: ", chinfo[:mixed_P])
    println("    mixed_S: ", chinfo[:mixed_S])
end
println()

println("uubar_to_ssbar:")
println("  类型: ", info2[:type])
for (ch, chinfo) in info2[:channels]
    println("  $ch 通道:")
    println("    simple: ", chinfo[:simple])
    println("    mixed_P: ", chinfo[:mixed_P])
    println("    mixed_S: ", chinfo[:mixed_S])
end
println()

# 对比散射振幅
println("="^80)
println("散射振幅 M² 对比")
println("="^80)
println()

# 在几个(s,t)点计算M²
s_threshold_ssbar = (masses.s + masses.s)^2
s_threshold_uubar = (masses.u + masses.u)^2

println("阈值:")
@printf("  s_threshold(ssbar) = %.6f fm⁻² = %.2f MeV²\n", s_threshold_ssbar, s_threshold_ssbar * ħc_MeV_fm^2)
@printf("  s_threshold(uubar) = %.6f fm⁻² = %.2f MeV²\n", s_threshold_uubar, s_threshold_uubar * ħc_MeV_fm^2)
println()

# 测试点：在ssbar阈值附近
s_test = s_threshold_ssbar * 1.5
t_values = [-0.5, -0.2, -0.1, 0.0]  # 不同的t值

println("在 s = $(round(s_test * ħc_MeV_fm^2, digits=2)) MeV² 处:")
println("-"^80)
@printf("%-15s %20s %20s %15s\n", "t (MeV²)", "M²(ssbar→uubar)", "M²(uubar→ssbar)", "比值")
println("-"^80)

for t in t_values
    M2_ssbar_uubar = scattering_amplitude_squared(proc1, s_test, t, quark_params, thermo_params, K_coeffs)
    M2_uubar_ssbar = scattering_amplitude_squared(proc2, s_test, t, quark_params, thermo_params, K_coeffs)
    ratio = M2_ssbar_uubar / M2_uubar_ssbar
    @printf("%15.2f %20.6e %20.6e %15.3f\n", t * ħc_MeV_fm^2, M2_ssbar_uubar, M2_uubar_ssbar, ratio)
end
println()

# 对比截面
println("="^80)
println("总截面 σ(s) 对比")
println("="^80)
println()

s_values = [s_threshold_ssbar * 1.1, s_threshold_ssbar * 1.5, s_threshold_ssbar * 2.0, s_threshold_ssbar * 5.0]

println("-"^80)
@printf("%-20s %20s %20s %15s\n", "s (MeV²)", "σ(ssbar→uubar)", "σ(uubar→ssbar)", "比值")
println("-"^80)

for s in s_values
    σ_ssbar_uubar = total_cross_section(proc1, s, quark_params, thermo_params, K_coeffs)
    σ_uubar_ssbar = total_cross_section(proc2, s, quark_params, thermo_params, K_coeffs)
    ratio = σ_ssbar_uubar / σ_uubar_ssbar
    @printf("%20.2f %20.6e %20.6e %15.3f\n", s * ħc_MeV_fm^2, σ_ssbar_uubar, σ_uubar_ssbar, ratio)
end
println()

# 检查对称性
println("="^80)
println("对称性检查")
println("="^80)
println()

println("理论上，对于 A+B → C+D 和 C+D → A+B:")
println("  如果质量相同，则 M²(s,t) 应该相同")
println("  但这里质量不同：m_s ≠ m_u")
println()

println("Mandelstam变量关系:")
println("  对于 ssbar→uubar: s + t + u = 2m_s² + 2m_u²")
println("  对于 uubar→ssbar: s + t + u = 2m_u² + 2m_s²")
println("  所以 s 和 t 的物理意义相同")
println()

# 计算在相同(s,t)点的M²
s_test = s_threshold_ssbar * 2.0
t_test = -0.1

println("在相同的 (s,t) 点:")
@printf("  s = %.2f MeV²\n", s_test * ħc_MeV_fm^2)
@printf("  t = %.2f MeV²\n", t_test * ħc_MeV_fm^2)
println()

M2_1 = scattering_amplitude_squared(proc1, s_test, t_test, quark_params, thermo_params, K_coeffs)
M2_2 = scattering_amplitude_squared(proc2, s_test, t_test, quark_params, thermo_params, K_coeffs)

@printf("  M²(ssbar→uubar) = %.6e\n", M2_1)
@printf("  M²(uubar→ssbar) = %.6e\n", M2_2)
@printf("  比值 = %.3f\n", M2_1 / M2_2)
println()

if M2_1 / M2_2 > 2.0
    println("⚠️ 散射振幅比值异常大！")
    println("   这说明 M² 的计算可能有问题")
    println()
    println("可能的原因:")
    println("  1. 介子通道定义不对称")
    println("  2. 传播子计算有bug")
    println("  3. 耦合常数使用错误")
    println("  4. Mandelstam变量计算错误")
else
    println("✓ 散射振幅比值正常")
end

println()
println("="^80)
println("结论")
println("="^80)
println()

println("需要检查:")
println("  1. 为什么 M²(ssbar→uubar) 比 M²(uubar→ssbar) 大这么多？")
println("  2. 介子通道是否对称？")
println("  3. 是否有额外的对称因子？")
println()

println("="^80)
