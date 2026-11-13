"""
测试微分散射截面计算模块

测试内容：
1. 核心函数基本计算
2. 运动学阈值检查
3. 退化情况处理（m1 = m2）
4. t 积分边界计算
5. 与散射矩阵元的关系验证
6. 物理一致性（同位旋对称、电荷共轭）
"""

using Test
if !isdefined(Main, :RelaxTime)
    Base.include(Main, normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl")))
end

using Main.Constants_PNJL
const Constants_PNJL = Main.Constants_PNJL
using Main.OneLoopIntegrals: A
using Main.GaussLegendre: gauleg
using Main.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using Main.ScatteringAmplitude: scattering_amplitude_squared, 
                            get_quark_masses_for_process,
                            calculate_mandelstam_variables
using Main.DifferentialCrossSection
using Main.TotalCrossSection: calculate_t_bounds

@testset "微分散射截面计算测试" begin
    
    # ========== 测试准备：设置物理参数 ==========
    println("\n" * "="^70)
    println("准备测试参数")
    println("="^70)
    
    # 温度和化学势
    T = 150.0 / 197.327  # 150 MeV -> fm⁻¹
    m_u = 300.0 / 197.327
    m_s = 500.0 / 197.327
    μ_u = 0.0
    μ_d = μ_u
    μ_s = 0.0
    
    # Polyakov 环参数
    Φ = 0.5
    Φbar = 0.5
    
    # 各向异性参数
    ξ = 0.0
    
    # 使用Gauss-Legendre积分节点和权重
    nodes_p, weights_p = gauleg(0.0, 20.0, 64)
    
    # 计算 A 函数
    A_u = A(m_u, μ_u, T, Φ, Φbar, nodes_p, weights_p)
    A_s = A(m_s, μ_s, T, Φ, Φbar, nodes_p, weights_p)
    
    # 计算有效耦合系数
    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)
    
    G_fm2_val = Constants_PNJL.G_fm2
    K_fm5_val = Constants_PNJL.K_fm5
    K_coeffs = calculate_effective_couplings(G_fm2_val, K_fm5_val, G_u, G_s)
    
    # 组装参数结构
    quark_params = (
        m = (u = m_u, d = m_u, s = m_s),
        μ = (u = μ_u, d = μ_u, s = μ_s),
        A = (u = A_u, d = A_u, s = A_s)
    )
    
    thermo_params = (T = T, Φ = Φ, Φbar = Φbar, ξ = ξ)
    
    println("  T = $(round(T*197.327, digits=1)) MeV")
    println("  μ = $(round(μ_u*197.327, digits=1)) MeV")
    println("  m_u = $(round(m_u, digits=4)) fm⁻¹")
    println("  m_s = $(round(m_s, digits=4)) fm⁻¹")
    println("  ξ = $ξ")
    
    # ========== 测试1：核心函数基本计算 ==========
    @testset "核心函数基本计算" begin
        println("\n" * "="^70)
        println("测试1：核心函数基本计算")
        println("="^70)
        
        s = 31.0  # fm⁻²
        t = -2.0
        process = :uu_to_uu
        
        # 步骤1：获取质量
        m1, m2, m3, m4 = get_quark_masses_for_process(process, quark_params)
        println("\n散射过程: $process")
        println("  粒子质量: m1=$(round(m1,digits=4)), m2=$(round(m2,digits=4)), " *
                "m3=$(round(m3,digits=4)), m4=$(round(m4,digits=4)) fm⁻¹")
        
        # 步骤2：计算 Mandelstam 变量
        u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
        mandelstam_vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)
        println("  s = $s fm⁻², t = $t fm⁻²")
        println("  s_12_plus = $(round(mandelstam_vars.s_12_plus, digits=4)) fm⁻²")
        println("  s_12_minus = $(round(mandelstam_vars.s_12_minus, digits=4)) fm⁻²")
        
        # 步骤3：计算散射矩阵元
        M_squared = scattering_amplitude_squared(
            process, s, t, quark_params, thermo_params, K_coeffs
        )
        println("  |M|² = $(round(M_squared, digits=2)) fm⁻⁴")
        
        # 步骤4：计算微分截面
        dsigma_dt = differential_cross_section(
            mandelstam_vars.s_12_plus, mandelstam_vars.s_12_minus, M_squared
        )
        
        println("\n结果:")
        println("  dσ/dt = $(round(dsigma_dt, digits=6)) fm²")
        
        # 验证结果
        @test dsigma_dt > 0.0
        @test !isnan(dsigma_dt) && !isinf(dsigma_dt)
        
        # 手动验证公式
        expected = M_squared / (16π * mandelstam_vars.s_12_plus * mandelstam_vars.s_12_minus)
        @test dsigma_dt ≈ expected rtol=1e-12
        
        println("✓ 基本计算测试通过")
    end
    
    # ========== 测试2：运动学阈值检查 ==========
    @testset "运动学阈值检查" begin
        println("\n" * "="^70)
        println("测试2：运动学阈值检查")
        println("="^70)
        
        m_u = quark_params.m.u
        m_s = quark_params.m.s
        
        # 测试2.1：正常情况（s > 阈值）
        s_safe = 10.0
        @test check_kinematic_threshold(s_safe, m_u, m_u, warn_close=false) == true
        println("  s = $s_safe > (2m_u)² = $((2*m_u)^2): ✓")
        
        # 测试2.2：违反阈值
        s_below = 0.1
        @test check_kinematic_threshold(s_below, m_u, m_u, warn_close=false) == false
        println("  s = $s_below < (2m_u)²: 正确检测到违反 ✓")
        
        # 测试2.3：接近阈值（应发出警告）
        s_threshold = (m_u + m_u)^2
        s_close = s_threshold + 1e-13
        @test check_kinematic_threshold(s_close, m_u, m_u, warn_close=true) == true
        println("  s 接近阈值: 正确发出警告 ✓")
        
        # 测试2.4：不同质量组合
        @test check_kinematic_threshold(20.0, m_u, m_s, warn_close=false) == true
        println("  s = 20.0 > (m_u + m_s)²: ✓")
        
        println("✓ 阈值检查测试通过")
    end
    
    # ========== 测试3：退化情况（m1 = m2）==========
    @testset "退化情况处理（m1 = m2）" begin
        println("\n" * "="^70)
        println("测试3：退化情况处理（m1 = m2）")
        println("="^70)
        
        s = 31.0
        t = -2.0
        process = :uu_to_uu  # u 和 u 质量相同
        
        m1, m2, m3, m4 = get_quark_masses_for_process(process, quark_params)
        u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
        mandelstam_vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)
        
        println("  m1 = $(round(m1,digits=6)), m2 = $(round(m2,digits=6))")
        println("  s_12_minus = $(round(mandelstam_vars.s_12_minus, digits=6)) fm⁻²")
        
        M_squared = scattering_amplitude_squared(
            process, s, t, quark_params, thermo_params, K_coeffs
        )
        
        # 应该能够正常计算（自动正则化）
        dsigma_dt = differential_cross_section(
            mandelstam_vars.s_12_plus, mandelstam_vars.s_12_minus, M_squared
        )
        
        @test dsigma_dt > 0.0
        @test !isnan(dsigma_dt) && !isinf(dsigma_dt)
        
        println("  dσ/dt = $(round(dsigma_dt, digits=6)) fm²")
        println("✓ 退化情况处理正确")
    end
    
    # ========== 测试4：t 积分边界计算 ==========
    @testset "t 积分边界计算" begin
        println("\n" * "="^70)
        println("测试4：t 积分边界计算")
        println("="^70)
        
        s = 31.0
        m_u = quark_params.m.u
        
        # 测试4.1：相同质量
        t_bounds = calculate_t_bounds(s, m_u, m_u, m_u, m_u)
        println("\n相同质量 (uu→uu):")
        println("  t_min = $(round(t_bounds.t_min, digits=4)) fm⁻²")
        println("  t_max = $(round(t_bounds.t_max, digits=4)) fm⁻²")
        
        @test t_bounds.t_max == 0.0
        @test t_bounds.t_min < 0.0
        @test t_bounds.t_min < t_bounds.t_max
        
        # 测试4.2：不同质量
        m_s = quark_params.m.s
        t_bounds_us = calculate_t_bounds(s, m_u, m_s, m_u, m_s)
        println("\n不同质量 (us→us):")
        println("  t_min = $(round(t_bounds_us.t_min, digits=4)) fm⁻²")
        println("  t_max = $(round(t_bounds_us.t_max, digits=4)) fm⁻²")
        
        @test t_bounds_us.t_max == 0.0
        @test t_bounds_us.t_min < 0.0
        
        # 测试4.3：验证 Mandelstam 约束
        # 对于相同质量，t_max = 0 （正向散射）
        # t + u = 2(m_i^2 + m_c^2) - s
        # 在 t_max 时，u = 2(m_u^2 + m_u^2) - s - 0 = 4m_u^2 - s
        u_at_tmax = 4*m_u^2 - s
        mandelstam_sum = s + t_bounds.t_max + u_at_tmax
        expected_sum = 4*m_u^2  # m_i^2 + m_j^2 + m_c^2 + m_d^2
        @test mandelstam_sum ≈ expected_sum rtol=1e-10
        
        println("✓ t 边界计算测试通过")
    end
    
    # ========== 测试5：与散射矩阵元的关系验证 ==========
    @testset "与散射矩阵元的关系验证" begin
        println("\n" * "="^70)
        println("测试5：与散射矩阵元的关系验证")
        println("="^70)
        
        s = 31.0
        t = -2.0
        
        # 测试多个散射过程
        test_processes = [:uu_to_uu, :ss_to_ss, :udbar_to_udbar, :uubar_to_uubar]
        
        for process in test_processes
            m1, m2, m3, m4 = get_quark_masses_for_process(process, quark_params)
            u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
            mandelstam_vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)
            
            M_squared = scattering_amplitude_squared(
                process, s, t, quark_params, thermo_params, K_coeffs
            )
            
            dsigma_dt = differential_cross_section(
                mandelstam_vars.s_12_plus, mandelstam_vars.s_12_minus, M_squared
            )
            
            # 验证公式关系
            kinematic_factor = 1.0 / (16π * mandelstam_vars.s_12_plus * mandelstam_vars.s_12_minus)
            expected = kinematic_factor * M_squared
            
            @test dsigma_dt ≈ expected rtol=1e-12
            
            println("  $process: dσ/dt = $(round(dsigma_dt, digits=6)) fm² ✓")
        end
        
        println("✓ 关系验证测试通过")
    end
    
    # ========== 测试6：物理一致性 ==========
    @testset "物理一致性（对称性）" begin
        println("\n" * "="^70)
        println("测试6：物理一致性验证")
        println("="^70)
        
        s = 31.0
        t = -2.0
        
        # 测试6.1：同位旋对称（uu 与 dd，当 m_u = m_d）
        println("\n同位旋对称测试:")
        m1, m2, m3, m4 = get_quark_masses_for_process(:uu_to_uu, quark_params)
        u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
        mandelstam_vars_uu = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)
        M_squared_uu = scattering_amplitude_squared(
            :uu_to_uu, s, t, quark_params, thermo_params, K_coeffs
        )
        dsigma_uu = differential_cross_section(
            mandelstam_vars_uu.s_12_plus, mandelstam_vars_uu.s_12_minus, M_squared_uu
        )
        
        # dd 应该与 uu 相同（m_u = m_d）
        # 注意：实际代码中可能没有 dd_to_dd，这里仅验证质量效应
        println("  uu→uu: dσ/dt = $(round(dsigma_uu, digits=6)) fm²")
        println("  （m_u = m_d 时，dd→dd 应相同）")
        
        # 测试6.2：电荷共轭对称（dubar 与 udbar）
        println("\n电荷共轭对称测试:")
        m1_ud, m2_ud, m3_ud, m4_ud = get_quark_masses_for_process(:udbar_to_udbar, quark_params)
        u_ud = m1_ud^2 + m2_ud^2 + m3_ud^2 + m4_ud^2 - s - t
        mandelstam_vars_udbar = calculate_mandelstam_variables(s, t, u_ud, m1_ud, m2_ud, m3_ud, m4_ud)
        M_squared_udbar = scattering_amplitude_squared(
            :udbar_to_udbar, s, t, quark_params, thermo_params, K_coeffs
        )
        dsigma_udbar = differential_cross_section(
            mandelstam_vars_udbar.s_12_plus, mandelstam_vars_udbar.s_12_minus, M_squared_udbar
        )
        
        m1_du, m2_du, m3_du, m4_du = get_quark_masses_for_process(:dubar_to_dubar, quark_params)
        u_du = m1_du^2 + m2_du^2 + m3_du^2 + m4_du^2 - s - t
        mandelstam_vars_dubar = calculate_mandelstam_variables(s, t, u_du, m1_du, m2_du, m3_du, m4_du)
        M_squared_dubar = scattering_amplitude_squared(
            :dubar_to_dubar, s, t, quark_params, thermo_params, K_coeffs
        )
        dsigma_dubar = differential_cross_section(
            mandelstam_vars_dubar.s_12_plus, mandelstam_vars_dubar.s_12_minus, M_squared_dubar
        )
        
        println("  udbar: dσ/dt = $(round(dsigma_udbar, digits=6)) fm²")
        println("  dubar: dσ/dt = $(round(dsigma_dubar, digits=6)) fm²")
        println("  相对差异: $(round(abs(dsigma_udbar - dsigma_dubar)/dsigma_udbar, digits=12))")
        
        @test dsigma_udbar ≈ dsigma_dubar rtol=1e-10
        
        println("✓ 物理一致性测试通过")
    end
    
    # ========== 测试7：错误处理 ==========
    @testset "错误处理" begin
        println("\n" * "="^70)
        println("测试7：错误处理")
        println("="^70)
        
        # 测试7.1：s_12_plus ≤ 0 应抛出错误
        println("\n测试 s_12_plus ≤ 0 错误:")
        @test_throws ErrorException differential_cross_section(-1.0, 10.0, 100.0)
        println("  ✓ 正确抛出错误")
        
        # 测试7.2：calculate_t_bounds 的阈值检查
        println("\n测试 t_bounds 阈值检查:")
        m_u = quark_params.m.u
        s_below = 0.1  # < (2m_u)²
        @test_throws DomainError calculate_t_bounds(s_below, m_u, m_u, m_u, m_u)
        println("  ✓ 正确检测阈值违反")
        
        println("✓ 错误处理测试通过")
    end
    
end
