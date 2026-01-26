#!/usr/bin/env julia
"""
步骤 2: 测试插值方法和积分节点数的影响
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
println("步骤 2: 测试收敛性")
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
@printf("  T = %.4f fm⁻¹\n", T)
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
println("测试 A: σ(s) 缓存点数的影响")
println("="^80)

p_grid_15, p_w_15 = gauleg(0.0, 15.0, 64)

println("\n| N (缓存点数) | w(ssbar→uubar) (fm⁻¹) | 相对变化 |")
println("|-------------|----------------------|----------|")

N_values = [64, 128, 256, 512]
w_results = Float64[]

for (idx, N) in enumerate(N_values)
    cache = RelaxationTime.AverageScatteringRate.build_w0cdf_pchip_cache(
        :ssbar_to_uubar,
        quark_params_with_A,
        thermo_params,
        K_coeffs;
        N = N,
        design_p_nodes = 14,
        design_angle_nodes = 4,
        design_phi_nodes = 8,
        p_cutoff = Λ_inv_fm,
        n_sigma_points = 32
    )
    
    w = RelaxationTime.AverageScatteringRate.average_scattering_rate(
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
    
    push!(w_results, w)
    
    if idx == 1
        @printf("| %4d | %.6e | baseline |\n", N, w)
    else
        rel_change = abs(w - w_results[idx-1]) / w_results[idx-1] * 100
        @printf("| %4d | %.6e | %.3f%% |\n", N, w, rel_change)
    end
end

println("\n" * "="^80)
println("测试 B: 动量积分节点数的影响")
println("="^80)

# 使用固定的缓存
cache_fixed = RelaxationTime.AverageScatteringRate.build_w0cdf_pchip_cache(
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

println("\n| p_nodes | w(ssbar→uubar) (fm⁻¹) | 相对变化 |")
println("|---------|----------------------|----------|")

p_nodes_values = [32, 64, 128, 256]
w_results_p = Float64[]

for (idx, p_nodes) in enumerate(p_nodes_values)
    p_grid, p_w = gauleg(0.0, 15.0, p_nodes)
    
    w = RelaxationTime.AverageScatteringRate.average_scattering_rate(
        :ssbar_to_uubar,
        quark_params_with_A,
        thermo_params,
        K_coeffs;
        p_grid = p_grid,
        p_w = p_w,
        angle_nodes = 64,
        phi_nodes = 8,
        cs_cache = cache_fixed,
        apply_s_domain_cut = true,
        sigma_cutoff = Λ_inv_fm
    )
    
    push!(w_results_p, w)
    
    if idx == 1
        @printf("| %4d | %.6e | baseline |\n", p_nodes, w)
    else
        rel_change = abs(w - w_results_p[idx-1]) / w_results_p[idx-1] * 100
        @printf("| %4d | %.6e | %.3f%% |\n", p_nodes, w, rel_change)
    end
end

println("\n" * "="^80)
println("测试 C: 角度积分节点数的影响")
println("="^80)

println("\n| angle_nodes | w(ssbar→uubar) (fm⁻¹) | 相对变化 |")
println("|-------------|----------------------|----------|")

angle_nodes_values = [32, 64, 128, 256]
w_results_angle = Float64[]

for (idx, angle_nodes) in enumerate(angle_nodes_values)
    w = RelaxationTime.AverageScatteringRate.average_scattering_rate(
        :ssbar_to_uubar,
        quark_params_with_A,
        thermo_params,
        K_coeffs;
        p_grid = p_grid_15,
        p_w = p_w_15,
        angle_nodes = angle_nodes,
        phi_nodes = 8,
        cs_cache = cache_fixed,
        apply_s_domain_cut = true,
        sigma_cutoff = Λ_inv_fm
    )
    
    push!(w_results_angle, w)
    
    if idx == 1
        @printf("| %4d | %.6e | baseline |\n", angle_nodes, w)
    else
        rel_change = abs(w - w_results_angle[idx-1]) / w_results_angle[idx-1] * 100
        @printf("| %4d | %.6e | %.3f%% |\n", angle_nodes, w, rel_change)
    end
end

println("\n" * "="^80)
println("测试 D: 方位角积分节点数的影响")
println("="^80)

println("\n| phi_nodes | w(ssbar→uubar) (fm⁻¹) | 相对变化 |")
println("|-----------|----------------------|----------|")

phi_nodes_values = [4, 8, 16, 32]
w_results_phi = Float64[]

for (idx, phi_nodes) in enumerate(phi_nodes_values)
    w = RelaxationTime.AverageScatteringRate.average_scattering_rate(
        :ssbar_to_uubar,
        quark_params_with_A,
        thermo_params,
        K_coeffs;
        p_grid = p_grid_15,
        p_w = p_w_15,
        angle_nodes = 64,
        phi_nodes = phi_nodes,
        cs_cache = cache_fixed,
        apply_s_domain_cut = true,
        sigma_cutoff = Λ_inv_fm
    )
    
    push!(w_results_phi, w)
    
    if idx == 1
        @printf("| %4d | %.6e | baseline |\n", phi_nodes, w)
    else
        rel_change = abs(w - w_results_phi[idx-1]) / w_results_phi[idx-1] * 100
        @printf("| %4d | %.6e | %.3f%% |\n", phi_nodes, w, rel_change)
    end
end

println("\n" * "="^80)
println("结论")
println("="^80)

println("\n如果:")
println("  - 测试 A: 增加缓存点数后差异 < 1% → 插值不是问题")
println("  - 测试 B: 增加动量节点后差异 < 1% → 动量积分精度足够")
println("  - 测试 C: 增加角度节点后差异 > 5% → 角度积分是关键!")
println("  - 测试 D: 增加方位角节点后差异 > 5% → 方位角积分是关键!")

println("\n与 Fortran 对比:")
w_fortran = 6.074e-2
w_julia_best = 0.0  # 需要从上面的测试中获取最好的结果
@printf("  Fortran: %.6e fm⁻¹\n", w_fortran)
println("  Julia (最佳设置): 见上面的测试结果")

println("\n测试完成!")
