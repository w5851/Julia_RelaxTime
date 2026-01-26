#!/usr/bin/env julia
"""
尝试精确匹配 Fortran 的积分设置
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
println("精确匹配 Fortran 积分设置")
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
println("Fortran 积分设置")
println("="^80)
println("  npoint_rex = 64")
println("  yp: [0, 1] (Gauss-Legendre)")
println("  ypi: [0, π] (Gauss-Legendre)")
println("  scale_p = 10.0")
println("  半无穷积分: p = scale_p * t / (1-t)")

println("\n" * "="^80)
println("构建 σ(s) 缓存")
println("="^80)

# 使用有限动量网格构建缓存 (与 Fortran 一致)
cache = RelaxationTime.AverageScatteringRate.build_w0cdf_pchip_cache(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    N = 128,  # 增加点数以提高精度
    design_p_nodes = 14,
    design_angle_nodes = 4,
    design_phi_nodes = 8,
    p_cutoff = Λ_inv_fm,
    n_sigma_points = 32
)

@printf("  σ(s) 范围: [%.4f, %.4f] fm⁻²\n", minimum(cache.s_vals), maximum(cache.s_vals))

println("\n" * "="^80)
println("测试不同的角度积分节点数")
println("="^80)

# Fortran 使用 1D 角度积分 (angel ∈ [0, π])
# Julia 使用 2D 角度积分 (cosθ_i, cosθ_j) + 1D 方位角 (φ)
# 
# 关键问题: Fortran 的 angel 是什么?
# - 如果是两粒子动量夹角,则应该对应 Julia 的 cosΘ
# - cosΘ = cosθ_i cosθ_j + sinθ_i sinθ_j cos(φ)
#
# 对于各向同性分布,积分 cosθ_i 和 cosθ_j 应该给出相同的结果

println("\n测试 1: 使用 Fortran 的节点数 (64)")
println("-"^80)

w1 = RelaxationTime.AverageScatteringRate.average_scattering_rate(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    p_nodes = 64,
    angle_nodes = 64,  # 对应 Fortran 的 ypi
    phi_nodes = 1,     # 尝试只用 1 个节点
    scale = 10.0,
    cs_cache = cache,
    apply_s_domain_cut = true,
    sigma_cutoff = Λ_inv_fm
)

@printf("  w(ssbar→uubar) = %.6e fm⁻¹\n", w1)

println("\n测试 2: 增加 phi 节点数")
println("-"^80)

w2 = RelaxationTime.AverageScatteringRate.average_scattering_rate(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    p_nodes = 64,
    angle_nodes = 64,
    phi_nodes = 8,
    scale = 10.0,
    cs_cache = cache,
    apply_s_domain_cut = true,
    sigma_cutoff = Λ_inv_fm
)

@printf("  w(ssbar→uubar) = %.6e fm⁻¹\n", w2)

println("\n测试 3: 减少角度节点数")
println("-"^80)

w3 = RelaxationTime.AverageScatteringRate.average_scattering_rate(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    p_nodes = 64,
    angle_nodes = 16,
    phi_nodes = 8,
    scale = 10.0,
    cs_cache = cache,
    apply_s_domain_cut = true,
    sigma_cutoff = Λ_inv_fm
)

@printf("  w(ssbar→uubar) = %.6e fm⁻¹\n", w3)

println("\n" * "="^80)
println("与 Fortran 对比")
println("="^80)

w_fortran = 6.124e-2  # fm⁻¹

@printf("\nFortran:  %.6e fm⁻¹\n", w_fortran)
@printf("Julia 1:  %.6e fm⁻¹ (p=64, angle=64, phi=1)\n", w1)
@printf("Julia 2:  %.6e fm⁻¹ (p=64, angle=64, phi=8)\n", w2)
@printf("Julia 3:  %.6e fm⁻¹ (p=64, angle=16, phi=8)\n", w3)

println("\n比值:")
@printf("  Julia 1 / Fortran = %.4f\n", w1 / w_fortran)
@printf("  Julia 2 / Fortran = %.4f\n", w2 / w_fortran)
@printf("  Julia 3 / Fortran = %.4f\n", w3 / w_fortran)

println("\n" * "="^80)
println("检查其他散射过程")
println("="^80)

# 测试 uubar→ssbar (应该一致)
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
    p_nodes = 64,
    angle_nodes = 64,
    phi_nodes = 8,
    scale = 10.0,
    cs_cache = cache_uubar_ssbar,
    apply_s_domain_cut = true,
    sigma_cutoff = Λ_inv_fm
)

w_uubar_ssbar_fortran = 4.724e-2  # fm⁻¹

println("\nuubar→ssbar:")
@printf("  Fortran: %.6e fm⁻¹\n", w_uubar_ssbar_fortran)
@printf("  Julia:   %.6e fm⁻¹\n", w_uubar_ssbar)
@printf("  比值:    %.4f\n", w_uubar_ssbar / w_uubar_ssbar_fortran)
@printf("  差异:    %.2f%%\n", abs(w_uubar_ssbar - w_uubar_ssbar_fortran) / w_uubar_ssbar_fortran * 100)

println("\n测试完成!")
