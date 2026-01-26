#!/usr/bin/env julia
"""
测试混合积分方案: 动量半无穷,但 s 截断

验证 Julia 使用 sigma_cutoff 参数后,是否与 Fortran 的混合方案一致
"""

using Printf

# 包含必要的模块
include("../../src/Constants_PNJL.jl")
using .Constants_PNJL

include("../../src/relaxtime/RelaxationTime.jl")
using .RelaxationTime

include("../../src/relaxtime/EffectiveCouplings.jl")
using .EffectiveCouplings

include("../../src/relaxtime/OneLoopIntegrals.jl")
using .OneLoopIntegrals

println("="^80)
println("测试混合积分方案: 动量半无穷 + s 截断")
println("="^80)

# 参数 (与 Fortran 一致)
const ħc = 0.1973269804  # GeV·fm
T_MeV = 300.0
μ_MeV = 0.67
T = T_MeV / (ħc * 1000.0)
μ = μ_MeV / (ħc * 1000.0)

Φ = 0.8404
m_u = 0.04051  # fm⁻¹
m_s = 1.0323   # fm⁻¹

println("\n参数:")
@printf("  T = %.4f fm⁻¹ (%.1f MeV)\n", T, T_MeV)
@printf("  μ = %.6f fm⁻¹ (%.2f MeV)\n", μ, μ_MeV)
@printf("  Φ = %.4f\n", Φ)
@printf("  m_u = %.5f fm⁻¹\n", m_u)
@printf("  m_s = %.5f fm⁻¹\n", m_s)
@printf("  Λ = %.4f fm⁻¹ (%.1f MeV)\n", Λ_inv_fm, Λ_inv_fm * ħc * 1000)

# 准备参数
quark_params = (
    m = (u = m_u, d = m_u, s = m_s),
    μ = (u = μ, d = μ, s = μ)
)

thermo_params = (
    T = T,
    Φ = Φ,
    Φbar = Φ,
    ξ = 0.0
)

# 计算 A 和 K
nodes_p, weights_p = RelaxationTime.AverageScatteringRate.gauleg(0.0, 20.0, 16)
A_u = OneLoopIntegrals.A(m_u, μ, T, Φ, Φ, nodes_p, weights_p)
A_s = OneLoopIntegrals.A(m_s, μ, T, Φ, Φ, nodes_p, weights_p)

G_u = EffectiveCouplings.calculate_G_from_A(A_u, m_u)
G_s = EffectiveCouplings.calculate_G_from_A(A_s, m_s)

K_coeffs = EffectiveCouplings.calculate_effective_couplings(
    Constants_PNJL.G_fm2,
    Constants_PNJL.K_fm5,
    G_u,
    G_s
)

quark_params_with_A = merge(quark_params, (A = (u = A_u, d = A_u, s = A_s),))

println("\n" * "="^80)
println("方案对比")
println("="^80)

# 方案 1: 默认 (半无穷动量,无 s 截断)
println("\n方案 1: 半无穷动量积分,无 s 截断 (Julia 默认)")
println("-"^80)

w_ssbar_uubar_default = RelaxationTime.AverageScatteringRate.average_scattering_rate(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    p_nodes = 20,
    angle_nodes = 4,
    phi_nodes = 8,
    scale = 10.0,
    sigma_cutoff = nothing  # 无截断
)

@printf("  w(ssbar→uubar) = %.6e fm⁻¹\n", w_ssbar_uubar_default)

# 方案 2: 使用 sigma_cutoff = Λ (混合方案)
println("\n方案 2: 半无穷动量积分,s 截断在 Λ (混合方案)")
println("-"^80)

w_ssbar_uubar_hybrid = RelaxationTime.AverageScatteringRate.average_scattering_rate(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    p_nodes = 20,
    angle_nodes = 4,
    phi_nodes = 8,
    scale = 10.0,
    sigma_cutoff = Λ_inv_fm  # 使用 Λ 截断
)

@printf("  w(ssbar→uubar) = %.6e fm⁻¹\n", w_ssbar_uubar_hybrid)

# 计算 s_up (Fortran 方式)
s_up_initial = (sqrt(m_s^2 + Λ_inv_fm^2) + sqrt(m_s^2 + Λ_inv_fm^2))^2
s_up_final = (sqrt(m_u^2 + Λ_inv_fm^2) + sqrt(m_u^2 + Λ_inv_fm^2))^2
s_up_fortran = min(s_up_initial, s_up_final)

@printf("\n  Fortran s_up = %.4f fm⁻² (min of initial=%.4f, final=%.4f)\n", 
        s_up_fortran, s_up_initial, s_up_final)

# 对比
println("\n" * "="^80)
println("结果对比")
println("="^80)

@printf("\n方案 1 (无截断):  %.6e fm⁻¹\n", w_ssbar_uubar_default)
@printf("方案 2 (Λ截断):   %.6e fm⁻¹\n", w_ssbar_uubar_hybrid)
@printf("比值 (方案1/方案2): %.4f\n", w_ssbar_uubar_default / w_ssbar_uubar_hybrid)
@printf("差异: %.2f%%\n", abs(w_ssbar_uubar_default - w_ssbar_uubar_hybrid) / w_ssbar_uubar_hybrid * 100)

println("\n" * "="^80)
println("与 Fortran 对比")
println("="^80)

# Fortran 结果 (从之前的输出)
w_fortran = 5.941e-2  # fm⁻¹

@printf("\nFortran (混合方案):  %.6e fm⁻¹\n", w_fortran)
@printf("Julia 方案 2:        %.6e fm⁻¹\n", w_ssbar_uubar_hybrid)
@printf("比值 (Julia/Fortran): %.4f\n", w_ssbar_uubar_hybrid / w_fortran)
@printf("差异: %.2f%%\n", abs(w_ssbar_uubar_hybrid - w_fortran) / w_fortran * 100)

if abs(w_ssbar_uubar_hybrid - w_fortran) / w_fortran < 0.02
    println("\n✅ 混合方案成功! Julia 与 Fortran 一致 (差异 < 2%)")
else
    println("\n⚠️  仍有差异,需要进一步调查")
end

println("\n" * "="^80)
println("物理解释")
println("="^80)

println("\n混合方案的物理意义:")
println("  1. 动量半无穷积分: 包含所有热分布的粒子")
println("  2. s 截断: 限制质心能量在模型有效范围内")
println("  3. 允许高动量的前向散射 (p > Λ, 小角度)")
println("  4. 排除高能量的大角度散射 (s > s_up)")

println("\n与纯动量截断的区别:")
println("  - 纯动量截断: p < Λ → 排除所有 p > Λ 的贡献")
println("  - 混合方案: s < s_up → 只排除高能量散射")
println("  - 混合方案包含更多物理贡献")

println("\n测试完成!")
