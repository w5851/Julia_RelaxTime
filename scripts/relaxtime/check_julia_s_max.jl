#!/usr/bin/env julia
"""
检查 Julia 在计算 ssbar→uubar 时使用的 s 范围
"""

using Printf

# 包含必要的模块
include("../../src/Constants_PNJL.jl")
using .Constants_PNJL

include("../../src/relaxtime/TotalCrossSection.jl")
using .TotalCrossSection

println("="^80)
println("检查 Julia 的 s_max 计算")
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
@printf("  Λ = %.4f fm⁻¹ (%.1f MeV)\n", Λ_inv_fm, Λ_inv_fm * ħc * 1000)

# 计算 s_threshold
s_th_ssbar_uubar = (2 * m_s)^2
s_th_uubar_ssbar = (2 * m_u)^2

println("\n阈值:")
@printf("  s_th(ssbar→uubar) = %.4f fm⁻² (%.1f MeV²)\n", s_th_ssbar_uubar, s_th_ssbar_uubar * (ħc * 1000)^2)
@printf("  s_th(uubar→ssbar) = %.6f fm⁻² (%.1f MeV²)\n", s_th_uubar_ssbar, s_th_uubar_ssbar * (ħc * 1000)^2)

# 计算 s_max (Fortran 方式)
s_max_initial = (sqrt(m_s^2 + Λ_inv_fm^2) + sqrt(m_s^2 + Λ_inv_fm^2))^2
s_max_final = (sqrt(m_u^2 + Λ_inv_fm^2) + sqrt(m_u^2 + Λ_inv_fm^2))^2
s_max_fortran = min(s_max_initial, s_max_final)

println("\ns_max 计算 (Fortran 方式):")
@printf("  初态 (s,sbar): (√(m_s²+Λ²) + √(m_s²+Λ²))² = %.4f fm⁻²\n", s_max_initial)
@printf("  末态 (u,ubar): (√(m_u²+Λ²) + √(m_u²+Λ²))² = %.4f fm⁻²\n", s_max_final)
@printf("  s_max = min(初态, 末态) = %.4f fm⁻²\n", s_max_fortran)

# 检查 Julia 的实现
println("\n" * "="^80)
println("检查 Julia 的 TotalCrossSection 模块")
println("="^80)

# 查看 build_w0cdf_pchip_cache 的实现
println("\n从代码中,我们知道:")
println("  1. p_cutoff 默认为 sigma_cutoff")
println("  2. sigma_cutoff 默认为 Λ_inv_fm")
println("  3. 在 build_w0cdf_pchip_cache 中:")
println("     - 使用 p_cutoff 来限制动量范围")
println("     - 但 s 的范围是如何确定的?")

# 尝试直接调用看看
println("\n尝试构建 σ(s) 缓存...")

try
    include("../../src/relaxtime/RelaxationTime.jl")
    using .RelaxationTime
    
    include("../../src/relaxtime/EffectiveCouplings.jl")
    using .EffectiveCouplings
    
    include("../../src/relaxtime/OneLoopIntegrals.jl")
    using .OneLoopIntegrals
    
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
    
    # 构建 σ(s) 缓存
    println("\n构建 ssbar→uubar 的 σ(s) 缓存...")
    cache = RelaxationTime.AverageScatteringRate.build_w0cdf_pchip_cache(
        :ssbar_to_uubar,
        quark_params_with_A,
        thermo_params,
        K_coeffs;
        p_cutoff = Λ_inv_fm,
        n_sigma_points = 64
    )
    
    println("\n缓存信息:")
    @printf("  s 范围: [%.4f, %.4f] fm⁻²\n", minimum(cache.s_vals), maximum(cache.s_vals))
    @printf("  σ 范围: [%.6e, %.6e] fm²\n", minimum(cache.sigma_vals), maximum(cache.sigma_vals))
    @printf("  点数: %d\n", length(cache.s_vals))
    
    # 对比
    println("\n" * "="^80)
    println("对比")
    println("="^80)
    
    s_max_julia = maximum(cache.s_vals)
    @printf("\nFortran s_max: %.4f fm⁻²\n", s_max_fortran)
    @printf("Julia s_max:   %.4f fm⁻²\n", s_max_julia)
    @printf("差异: %.2f%%\n", abs(s_max_julia - s_max_fortran) / s_max_fortran * 100)
    
    if abs(s_max_julia - s_max_fortran) / s_max_fortran < 0.01
        println("\n✅ s_max 一致!")
    else
        println("\n❌ s_max 不一致!")
        println("这可能是差异的来源!")
    end
    
    # 检查在 s > s_max_fortran 的 σ 值
    s_beyond = cache.s_vals[cache.s_vals .> s_max_fortran]
    if !isempty(s_beyond)
        println("\n⚠️  Julia 包含了 s > s_max_fortran 的点:")
        @printf("  数量: %d\n", length(s_beyond))
        @printf("  范围: [%.4f, %.4f] fm⁻²\n", minimum(s_beyond), maximum(s_beyond))
        
        # 这些点的 σ 值
        idx_beyond = findall(cache.s_vals .> s_max_fortran)
        sigma_beyond = cache.sigma_vals[idx_beyond]
        @printf("  σ 范围: [%.6e, %.6e] fm²\n", minimum(sigma_beyond), maximum(sigma_beyond))
        @printf("  平均 σ: %.6e fm²\n", sum(sigma_beyond) / length(sigma_beyond))
        
        println("\n这些额外的贡献可能解释了差异!")
    else
        println("\n✅ Julia 没有超出 Fortran 的 s 范围")
    end
    
catch e
    println("\n❌ 计算失败:")
    println(e)
    if isa(e, ErrorException)
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
    end
end

println("\n分析完成!")
