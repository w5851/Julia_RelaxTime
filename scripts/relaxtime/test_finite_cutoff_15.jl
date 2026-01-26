#!/usr/bin/env julia
"""
测试使用有限截断 p_max = 15 fm⁻¹ 模拟半无穷积分
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
println("测试有限截断 p_max = 15 fm⁻¹")
println("="^80)

# 参数
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
@printf("  Λ = %.4f fm⁻¹\n", Λ_inv_fm)

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
println("构建 σ(s) 缓存 (使用 Λ 截断)")
println("="^80)

cache = RelaxationTime.AverageScatteringRate.build_w0cdf_pchip_cache(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    N = 128,
    design_p_nodes = 14,
    design_angle_nodes = 4,
    design_phi_nodes = 8,
    p_cutoff = Λ_inv_fm,
    n_sigma_points = 32
)

@printf("  σ(s) 范围: [%.4f, %.4f] fm⁻²\n", minimum(cache.s_vals), maximum(cache.s_vals))

println("\n" * "="^80)
println("测试不同的动量截断")
println("="^80)

p_max_values = [Λ_inv_fm, 10.0, 15.0, 20.0]

println("\n| p_max (fm⁻¹) | w(ssbar→uubar) (fm⁻¹) | 相对于 p=15 |")
println("|--------------|----------------------|-------------|")

w_15 = 0.0
for (idx, p_max) in enumerate(p_max_values)
    p_grid, p_w = RelaxationTime.AverageScatteringRate.gauleg(0.0, p_max, 64)
    
    w = RelaxationTime.AverageScatteringRate.average_scattering_rate(
        :ssbar_to_uubar,
        quark_params_with_A,
        thermo_params,
        K_coeffs;
        p_grid = p_grid,
        p_w = p_w,
        angle_nodes = 64,
        phi_nodes = 8,
        cs_cache = cache,
        apply_s_domain_cut = true,
        sigma_cutoff = Λ_inv_fm
    )
    
    if p_max == 15.0
        global w_15 = w
        @printf("| %.1f | %.6e | baseline |\n", p_max, w)
    else
        rel_diff = abs(w - w_15) / w_15 * 100
        @printf("| %.1f | %.6e | %.2f%% |\n", p_max, w, rel_diff)
    end
end

println("\n" * "="^80)
println("与 Fortran 对比 (p_max = 15 fm⁻¹)")
println("="^80)

p_grid_15, p_w_15 = RelaxationTime.AverageScatteringRate.gauleg(0.0, 15.0, 64)
w_julia = RelaxationTime.AverageScatteringRate.average_scattering_rate(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    p_grid = p_grid_15,
    p_w = p_w_15,
    angle_nodes = 64,
    phi_nodes = 8,
    cs_cache = cache,
    apply_s_domain_cut = true,
    sigma_cutoff = Λ_inv_fm
)

w_fortran = 6.074e-2  # fm⁻¹

@printf("\nFortran (p_max=15):  %.6e fm⁻¹\n", w_fortran)
@printf("Julia (p_max=15):    %.6e fm⁻¹\n", w_julia)
@printf("比值 (Julia/Fortran): %.4f\n", w_julia / w_fortran)
@printf("差异: %.2f%%\n", abs(w_julia - w_fortran) / w_fortran * 100)

if abs(w_julia - w_fortran) / w_fortran < 0.05
    println("\n✅ 成功! Julia 与 Fortran 一致 (差异 < 5%)")
elseif abs(w_julia - w_fortran) / w_fortran < 0.10
    println("\n⚠️  接近一致 (差异 5-10%)")
else
    println("\n❌ 仍有显著差异 (> 10%)")
    println("\n可能的原因:")
    println("  1. 角度积分方式不同 (Fortran: 1D angel, Julia: 2D cosθ + φ)")
    println("  2. 积分节点数不同")
    println("  3. σ(s) 插值方法不同")
end

println("\n" * "="^80)
println("检查其他散射过程")
println("="^80)

cache_uubar_ssbar = RelaxationTime.AverageScatteringRate.build_w0cdf_pchip_cache(
    :uubar_to_ssbar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    N = 128,
    design_p_nodes = 14,
    design_angle_nodes = 4,
    design_phi_nodes = 8,
    p_cutoff = Λ_inv_fm,
    n_sigma_points = 32
)

w_uubar_ssbar = RelaxationTime.AverageScatteringRate.average_scattering_rate(
    :uubar_to_ssbar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    p_grid = p_grid_15,
    p_w = p_w_15,
    angle_nodes = 64,
    phi_nodes = 8,
    cs_cache = cache_uubar_ssbar,
    apply_s_domain_cut = true,
    sigma_cutoff = Λ_inv_fm
)

w_uubar_ssbar_fortran = 4.692e-2  # fm⁻¹

println("\nuubar→ssbar:")
@printf("  Fortran: %.6e fm⁻¹\n", w_uubar_ssbar_fortran)
@printf("  Julia:   %.6e fm⁻¹\n", w_uubar_ssbar)
@printf("  比值:    %.4f\n", w_uubar_ssbar / w_uubar_ssbar_fortran)
@printf("  差异:    %.2f%%\n", abs(w_uubar_ssbar - w_uubar_ssbar_fortran) / w_uubar_ssbar_fortran * 100)

if abs(w_uubar_ssbar - w_uubar_ssbar_fortran) / w_uubar_ssbar_fortran < 0.05
    println("  ✅ 一致!")
end

println("\n测试完成!")
