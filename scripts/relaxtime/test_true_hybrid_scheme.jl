#!/usr/bin/env julia
"""
测试真正的混合积分方案:
- σ(s) 缓存: 使用有限动量网格 [0, Λ] 确定 s 范围
- 实际积分: 使用半无穷动量积分,但应用 s 截断

这与 Fortran 的实现一致
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
println("测试真正的混合积分方案")
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
println("真正的混合方案实现")
println("="^80)

# 步骤 1: 使用有限动量网格构建 σ(s) 缓存
println("\n步骤 1: 构建 σ(s) 缓存 (使用有限动量网格 [0, Λ])")
println("-"^80)

cache = RelaxationTime.AverageScatteringRate.build_w0cdf_pchip_cache(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    N = 64,
    design_p_nodes = 14,
    design_angle_nodes = 4,
    design_phi_nodes = 8,
    p_cutoff = Λ_inv_fm,  # 使用有限动量网格
    n_sigma_points = 32
)

s_min = minimum(cache.s_vals)
s_max = maximum(cache.s_vals)
@printf("  σ(s) 缓存范围: [%.4f, %.4f] fm⁻²\n", s_min, s_max)
@printf("  缓存点数: %d\n", length(cache.s_vals))

# 计算 Fortran 的 s_up
s_up_initial = (sqrt(m_s^2 + Λ_inv_fm^2) + sqrt(m_s^2 + Λ_inv_fm^2))^2
s_up_final = (sqrt(m_u^2 + Λ_inv_fm^2) + sqrt(m_u^2 + Λ_inv_fm^2))^2
s_up_fortran = min(s_up_initial, s_up_final)

@printf("\n  Fortran s_up = %.4f fm⁻² (min of initial=%.4f, final=%.4f)\n", 
        s_up_fortran, s_up_initial, s_up_final)
@printf("  Julia s_max = %.4f fm⁻²\n", s_max)
@printf("  差异: %.2f%%\n", abs(s_max - s_up_fortran) / s_up_fortran * 100)

# 步骤 2: 使用半无穷动量积分,但应用 s 截断
println("\n步骤 2: 计算散射率 (半无穷动量积分 + s 截断)")
println("-"^80)

w_hybrid = RelaxationTime.AverageScatteringRate.average_scattering_rate(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    p_nodes = 20,
    angle_nodes = 4,
    phi_nodes = 8,
    scale = 10.0,
    cs_cache = cache,  # 使用预构建的缓存
    apply_s_domain_cut = true,  # 应用 s 截断
    sigma_cutoff = Λ_inv_fm  # 用于确定 s_up
)

@printf("  w(ssbar→uubar) = %.6e fm⁻¹\n", w_hybrid)

# 对比
println("\n" * "="^80)
println("与 Fortran 对比")
println("="^80)

# Fortran 结果
w_fortran = 5.941e-2  # fm⁻¹

@printf("\nFortran (混合方案):  %.6e fm⁻¹\n", w_fortran)
@printf("Julia (真正混合):    %.6e fm⁻¹\n", w_hybrid)
@printf("比值 (Julia/Fortran): %.4f\n", w_hybrid / w_fortran)
@printf("差异: %.2f%%\n", abs(w_hybrid - w_fortran) / w_fortran * 100)

if abs(w_hybrid - w_fortran) / w_fortran < 0.05
    println("\n✅ 成功! Julia 与 Fortran 一致 (差异 < 5%)")
elseif abs(w_hybrid - w_fortran) / w_fortran < 0.10
    println("\n⚠️  接近一致,但仍有小差异 (5-10%)")
else
    println("\n❌ 仍有显著差异 (> 10%)")
    println("\n可能的原因:")
    println("  1. 积分节点数不同")
    println("  2. scale_p 参数不同")
    println("  3. 数值积分方法的细微差异")
    println("  4. σ(s) 插值方法不同")
end

println("\n" * "="^80)
println("方案说明")
println("="^80)

println("\n真正的混合方案:")
println("  1. σ(s) 缓存构建:")
println("     - 使用有限动量网格 [0, Λ]")
println("     - 确定 s 的有效范围 [s_min, s_max]")
println("     - 预计算 σ(s) 并插值")
println()
println("  2. 散射率积分:")
println("     - 使用半无穷动量积分 [0, ∞)")
println("     - 对每个 (p_i, p_j, Θ) 计算 s")
println("     - 如果 s > s_max,则 σ(s) = 0 (缓存外)")
println("     - 如果 s < s_max,则使用插值的 σ(s)")
println()
println("  3. 物理意义:")
println("     - 允许高动量粒子 (p > Λ)")
println("     - 但限制质心能量 (s < s_max)")
println("     - 高动量的前向散射被包含")
println("     - 高能量的大角度散射被排除")

println("\n测试完成!")
