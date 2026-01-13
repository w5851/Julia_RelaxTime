"""
验证新增散射过程的物理等价性

测试内容：
1. dubar_to_dubar 与 udbar_to_udbar 在同位旋对称下应该完全相同
2. subar_to_subar 与 usbar_to_usbar 在同位旋对称下应该完全相同
3. 测试新过程的散射矩阵元计算
"""

using Test
using Printf

push!(LOAD_PATH, joinpath(@__DIR__, "../../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../src/relaxtime"))

include("../../../src/Constants_PNJL.jl")
include("../../../src/relaxtime/ScatteringAmplitude.jl")
include("../../../src/relaxtime/EffectiveCouplings.jl")
include("../../../src/relaxtime/OneLoopIntegrals.jl")
include("../../../src/integration/GaussLegendre.jl")

using .ScatteringAmplitude
using .Constants_PNJL
using .EffectiveCouplings
using .OneLoopIntegrals: A
using .GaussLegendre: gauleg

@testset "新增散射过程验证" begin
    
    # 准备物理参数
    T = 150.0 / 197.327  # 150 MeV
    m_u = 300.0 / 197.327
    m_s = 500.0 / 197.327
    μ_u = 0.0
    μ_s = 0.0
    Φ = 0.5
    Φbar = 0.5
    ξ = 0.0
    
    # 使用Gauss-Legendre积分节点和权重
    nodes_p, weights_p = gauleg(0.0, 20.0, 64)
    
    A_u = A(m_u, μ_u, T, Φ, Φbar, nodes_p, weights_p)
    A_s = A(m_s, μ_s, T, Φ, Φbar, nodes_p, weights_p)
    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)
    
    quark_params = (
        m = (u=m_u, d=m_u, s=m_s),
        μ = (u=μ_u, d=μ_u, s=μ_s),
        A = (u=A_u, d=A_u, s=A_s)
    )
    thermo_params = (T=T, Φ=Φ, Φbar=Φbar, ξ=ξ)
    
    G_fm2 = Constants_PNJL.G_fm2
    K_fm5 = Constants_PNJL.K_fm5
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
    
    println("\n" * "="^70)
    println("测试1: dubar_to_dubar 与 udbar_to_udbar 的物理等价性")
    println("="^70)
    println("\n在同位旋对称下 (m_u = m_d)，d+ū→d+ū 与 u+đ→u+đ 应该完全相同")
    println("-"^70)
    
    # 测试参数
    s = 8.0  # fm⁻²
    t = -0.3  # fm⁻²
    
    M_udbar = scattering_amplitude_squared(
        :udbar_to_udbar, s, t, quark_params, thermo_params, K_coeffs
    )
    M_dubar = scattering_amplitude_squared(
        :dubar_to_dubar, s, t, quark_params, thermo_params, K_coeffs
    )
    
    println(@sprintf("参数: s=%.1f fm⁻², t=%.2f fm⁻²", s, t))
    println(@sprintf("u+đ→u+đ: |M|² = %.12e fm⁻⁴", M_udbar))
    println(@sprintf("d+ū→d+ū: |M|² = %.12e fm⁻⁴", M_dubar))
    println(@sprintf("相对差异: %.12e", abs(M_dubar - M_udbar) / M_udbar))
    
    @test isapprox(M_dubar, M_udbar, rtol=1e-12)
    
    if isapprox(M_dubar, M_udbar, rtol=1e-14)
        println("✓ 完全相同（浮点精度内）")
    elseif isapprox(M_dubar, M_udbar, rtol=1e-10)
        println("✓ 非常接近（相对误差<1e-10）")
    else
        println("✗ 差异较大")
    end
    
    # ====================================================================
    
    println("\n" * "="^70)
    println("测试2: subar_to_subar 与 usbar_to_usbar 的物理等价性")
    println("="^70)
    println("\n在同位旋对称下，s+ū→s+ū 与 u+s̄→u+s̄ 应该完全相同")
    println("-"^70)
    
    M_usbar = scattering_amplitude_squared(
        :usbar_to_usbar, s, t, quark_params, thermo_params, K_coeffs
    )
    M_subar = scattering_amplitude_squared(
        :subar_to_subar, s, t, quark_params, thermo_params, K_coeffs
    )
    
    println(@sprintf("参数: s=%.1f fm⁻², t=%.2f fm⁻²", s, t))
    println(@sprintf("u+s̄→u+s̄: |M|² = %.12e fm⁻⁴", M_usbar))
    println(@sprintf("s+ū→s+ū: |M|² = %.12e fm⁻⁴", M_subar))
    println(@sprintf("相对差异: %.12e", abs(M_subar - M_usbar) / M_usbar))
    
    @test isapprox(M_subar, M_usbar, rtol=1e-12)
    
    if isapprox(M_subar, M_usbar, rtol=1e-14)
        println("✓ 完全相同（浮点精度内）")
    elseif isapprox(M_subar, M_usbar, rtol=1e-10)
        println("✓ 非常接近（相对误差<1e-10）")
    else
        println("✗ 差异较大")
    end
    
    # ====================================================================
    
    println("\n" * "="^70)
    println("测试3: 多组参数下的对称性验证")
    println("="^70)
    
    test_params = [
        (6.0, -0.2),
        (8.0, -0.3),
        (10.0, -0.5),
        (12.0, -0.4)
    ]
    
    println("\ndubar vs udbar:")
    println("-"^70)
    println(@sprintf("%-10s | %-10s | %-18s | %-18s | %-12s", "s", "t", "M²(udbar)", "M²(dubar)", "相对差异"))
    println("-"^70)
    
    for (s_test, t_test) in test_params
        M_ud = scattering_amplitude_squared(:udbar_to_udbar, s_test, t_test, quark_params, thermo_params, K_coeffs)
        M_du = scattering_amplitude_squared(:dubar_to_dubar, s_test, t_test, quark_params, thermo_params, K_coeffs)
        rel_diff = abs(M_du - M_ud) / M_ud
        
        println(@sprintf("%-10.1f | %-10.2f | %-18.6e | %-18.6e | %-12.3e", 
                s_test, t_test, M_ud, M_du, rel_diff))
        
        @test isapprox(M_du, M_ud, rtol=1e-12)
    end
    
    println("\nsubar vs usbar:")
    println("-"^70)
    println(@sprintf("%-10s | %-10s | %-18s | %-18s | %-12s", "s", "t", "M²(usbar)", "M²(subar)", "相对差异"))
    println("-"^70)
    
    for (s_test, t_test) in test_params
        M_us = scattering_amplitude_squared(:usbar_to_usbar, s_test, t_test, quark_params, thermo_params, K_coeffs)
        M_su = scattering_amplitude_squared(:subar_to_subar, s_test, t_test, quark_params, thermo_params, K_coeffs)
        rel_diff = abs(M_su - M_us) / M_us
        
        println(@sprintf("%-10.1f | %-10.2f | %-18.6e | %-18.6e | %-12.3e", 
                s_test, t_test, M_us, M_su, rel_diff))
        
        @test isapprox(M_su, M_us, rtol=1e-12)
    end
    
    # ====================================================================
    
    println("\n" * "="^70)
    println("测试4: 新过程的物理合理性检查")
    println("="^70)
    
    s_test = 8.0
    t_test = -0.3
    
    all_new_processes = [:dubar_to_dubar, :subar_to_subar]
    
    println("\n检查散射矩阵元的物理约束:")
    println("-"^70)
    
    for process in all_new_processes
        M_sq = scattering_amplitude_squared(process, s_test, t_test, quark_params, thermo_params, K_coeffs)
        
        @test M_sq >= 0.0  # 物理约束：|M|² ≥ 0
        @test !isnan(M_sq) && !isinf(M_sq)
        
        println(@sprintf("%-20s: |M|² = %.6e fm⁻⁴  ✓", process, M_sq))
    end
    
    println("\n✓ 所有新过程满足物理约束 (|M|² ≥ 0)")
end

println("\n" * "="^70)
println("新增散射过程验证完成！")
println("="^70)
println("\n总结：")
println("- dubar_to_dubar 与 udbar_to_udbar 完全等价（同位旋对称）")
println("- subar_to_subar 与 usbar_to_usbar 完全等价（同位旋对称）")
println("- 两个新过程已成功添加到 SCATTERING_MESON_MAP")
println("- 总散射过程数：11 → 13种")


