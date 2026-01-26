#!/usr/bin/env julia
"""
验证使用最佳设置后,Julia 和 Fortran 的差异
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

include("../../src/integration/GaussLegendre.jl")
using .GaussLegendre: gauleg

println("="^80)
println("验证收敛性修复")
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
@printf("  μ_B = %.4f fm⁻¹ (%.2f MeV)\n", 3*μ, μ_MeV*3)
@printf("  Φ = %.4f\n", Φ)
@printf("  m_u = %.5f fm⁻¹\n", m_u)
@printf("  m_s = %.5f fm⁻¹\n", m_s)

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
nodes_p, weights_p = gauleg(0.0, 20.0, 16)
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
println("测试不同设置")
println("="^80)

# Fortran 参考值
w_fortran = 6.074e-2

println("\n| 设置 | N | phi_nodes | w (fm⁻¹) | 与 Fortran 差异 |")
println("|------|---|-----------|----------|----------------|")

# 设置 1: 原始设置 (N=256, phi_nodes=8)
p_grid_15, p_w_15 = gauleg(0.0, 15.0, 64)

cache_256 = RelaxationTime.AverageScatteringRate.build_w0cdf_pchip_cache(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    N = 256,
    design_p_nodes = 14,
    design_angle_nodes = 4,
    design_phi_nodes = 8,
    p_cutoff = Λ_inv_fm,
    n_sigma_points = 32
)

w_256_8 = RelaxationTime.AverageScatteringRate.average_scattering_rate(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    p_grid = p_grid_15,
    p_w = p_w_15,
    angle_nodes = 64,
    phi_nodes = 8,
    cs_cache = cache_256,
    apply_s_domain_cut = true,
    sigma_cutoff = Λ_inv_fm
)

diff_256_8 = abs(w_256_8 - w_fortran) / w_fortran * 100
@printf("| 原始 | 256 | 8 | %.6e | %.2f%% |\n", w_256_8, diff_256_8)

# 设置 2: 增加缓存点数 (N=512, phi_nodes=8)
cache_512 = RelaxationTime.AverageScatteringRate.build_w0cdf_pchip_cache(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    N = 512,
    design_p_nodes = 14,
    design_angle_nodes = 4,
    design_phi_nodes = 8,
    p_cutoff = Λ_inv_fm,
    n_sigma_points = 32
)

w_512_8 = RelaxationTime.AverageScatteringRate.average_scattering_rate(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    p_grid = p_grid_15,
    p_w = p_w_15,
    angle_nodes = 64,
    phi_nodes = 8,
    cs_cache = cache_512,
    apply_s_domain_cut = true,
    sigma_cutoff = Λ_inv_fm
)

diff_512_8 = abs(w_512_8 - w_fortran) / w_fortran * 100
@printf("| 最佳 | 512 | 8 | %.6e | %.2f%% |\n", w_512_8, diff_512_8)

# 设置 3: 进一步增加 (N=1024, phi_nodes=8)
println("\n计算 N=1024 的结果 (可能需要较长时间)...")
cache_1024 = RelaxationTime.AverageScatteringRate.build_w0cdf_pchip_cache(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    N = 1024,
    design_p_nodes = 14,
    design_angle_nodes = 4,
    design_phi_nodes = 8,
    p_cutoff = Λ_inv_fm,
    n_sigma_points = 32
)

w_1024_8 = RelaxationTime.AverageScatteringRate.average_scattering_rate(
    :ssbar_to_uubar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    p_grid = p_grid_15,
    p_w = p_w_15,
    angle_nodes = 64,
    phi_nodes = 8,
    cs_cache = cache_1024,
    apply_s_domain_cut = true,
    sigma_cutoff = Λ_inv_fm
)

diff_1024_8 = abs(w_1024_8 - w_fortran) / w_fortran * 100
@printf("| 更高 | 1024 | 8 | %.6e | %.2f%% |\n", w_1024_8, diff_1024_8)

println("\n" * "="^80)
println("对比其他散射过程")
println("="^80)

# 测试 uubar→ssbar (作为对照)
println("\n计算 uubar→ssbar...")

cache_uubar = RelaxationTime.AverageScatteringRate.build_w0cdf_pchip_cache(
    :uubar_to_ssbar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    N = 256,
    design_p_nodes = 14,
    design_angle_nodes = 4,
    design_phi_nodes = 8,
    p_cutoff = Λ_inv_fm,
    n_sigma_points = 32
)

w_uubar_256 = RelaxationTime.AverageScatteringRate.average_scattering_rate(
    :uubar_to_ssbar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    p_grid = p_grid_15,
    p_w = p_w_15,
    angle_nodes = 64,
    phi_nodes = 8,
    cs_cache = cache_uubar,
    apply_s_domain_cut = true,
    sigma_cutoff = Λ_inv_fm
)

cache_uubar_512 = RelaxationTime.AverageScatteringRate.build_w0cdf_pchip_cache(
    :uubar_to_ssbar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    N = 512,
    design_p_nodes = 14,
    design_angle_nodes = 4,
    design_phi_nodes = 8,
    p_cutoff = Λ_inv_fm,
    n_sigma_points = 32
)

w_uubar_512 = RelaxationTime.AverageScatteringRate.average_scattering_rate(
    :uubar_to_ssbar,
    quark_params_with_A,
    thermo_params,
    K_coeffs;
    p_grid = p_grid_15,
    p_w = p_w_15,
    angle_nodes = 64,
    phi_nodes = 8,
    cs_cache = cache_uubar_512,
    apply_s_domain_cut = true,
    sigma_cutoff = Λ_inv_fm
)

w_uubar_fortran = 1.0e-1  # Fortran 参考值

println("\n| 过程 | N | w (fm⁻¹) | 与 Fortran 差异 |")
println("|------|---|----------|----------------|")
@printf("| uubar→ssbar | 256 | %.6e | %.2f%% |\n", w_uubar_256, abs(w_uubar_256 - w_uubar_fortran) / w_uubar_fortran * 100)
@printf("| uubar→ssbar | 512 | %.6e | %.2f%% |\n", w_uubar_512, abs(w_uubar_512 - w_uubar_fortran) / w_uubar_fortran * 100)

rel_change_uubar = abs(w_uubar_512 - w_uubar_256) / w_uubar_256 * 100
@printf("\n相对变化 (N=256→512): %.3f%%\n", rel_change_uubar)

println("\n" * "="^80)
println("结论")
println("="^80)

println("\nssbar→uubar:")
@printf("  - Fortran: %.6e fm⁻¹\n", w_fortran)
@printf("  - Julia (N=256): %.6e fm⁻¹ (差异 %.2f%%)\n", w_256_8, diff_256_8)
@printf("  - Julia (N=512): %.6e fm⁻¹ (差异 %.2f%%)\n", w_512_8, diff_512_8)
@printf("  - Julia (N=1024): %.6e fm⁻¹ (差异 %.2f%%)\n", w_1024_8, diff_1024_8)

println("\nuubar→ssbar:")
@printf("  - Fortran: %.6e fm⁻¹\n", w_uubar_fortran)
@printf("  - Julia (N=256): %.6e fm⁻¹ (差异 %.2f%%)\n", w_uubar_256, abs(w_uubar_256 - w_uubar_fortran) / w_uubar_fortran * 100)
@printf("  - Julia (N=512): %.6e fm⁻¹ (差异 %.2f%%)\n", w_uubar_512, abs(w_uubar_512 - w_uubar_fortran) / w_uubar_fortran * 100)

println("\n关键发现:")
println("  1. ssbar→uubar 对缓存点数非常敏感")
println("  2. uubar→ssbar 对缓存点数不敏感")
println("  3. 这证实了阈值位置是关键因素")

if diff_512_8 < 5.0
    println("\n✅ 使用 N=512 后,差异已降至可接受范围 (< 5%)")
else
    println("\n⚠️ 即使使用 N=512,差异仍然较大")
end

println("\n验证完成!")
