"""
测试总传播子计算模块

测试内容：
1. 味因子查询功能
2. 散射过程解析
3. 味因子自动确定（t/u/s道）
4. 质量提取功能（新增）
5. 质心系动量计算（新增）
6. 通道分离接口（新增）
7. 所有11种散射过程（新增）
8. 错误处理
"""

using Test

push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../src/relaxtime"))

include("../src/relaxtime/TotalPropagator.jl")
include("../src/Constants_PNJL.jl")
include("../src/relaxtime/EffectiveCouplings.jl")
include("../src/relaxtime/OneLoopIntegrals.jl")
using .TotalPropagator
using .Constants_PNJL
using .EffectiveCouplings
using .OneLoopIntegrals: A
using .TotalPropagator: extract_quark_flavor, parse_scattering_process, get_flavor_factors_for_channel

@testset "总传播子计算测试" begin
    
    # 测试1：味提取功能
    @testset "味提取功能" begin
        @test extract_quark_flavor(:u) == :u
        @test extract_quark_flavor(:d) == :d
        @test extract_quark_flavor(:s) == :s
        @test extract_quark_flavor(:ubar) == :u
        @test extract_quark_flavor(:dbar) == :d
        @test extract_quark_flavor(:sbar) == :s
    end
    
    # 测试2：味因子查询（表5.3）
    @testset "味因子查询" begin
        # 对角元
        @test get_flavor_factor(:u, :u) == 1.0
        @test get_flavor_factor(:d, :d) == -1.0
        @test get_flavor_factor(:s, :s) == 2.0
        
        # 非对角元
        @test get_flavor_factor(:u, :d) ≈ √2
        @test get_flavor_factor(:d, :u) ≈ √2
        @test get_flavor_factor(:u, :s) ≈ √2
        @test get_flavor_factor(:s, :u) ≈ √2
        @test get_flavor_factor(:d, :s) ≈ √2
        @test get_flavor_factor(:s, :d) ≈ √2
        
        # 对称性
        @test get_flavor_factor(:u, :d) == get_flavor_factor(:d, :u)
        @test get_flavor_factor(:u, :s) == get_flavor_factor(:s, :u)
        
        # 反夸克视为同味
        @test get_flavor_factor(:ubar, :u) == 1.0
        @test get_flavor_factor(:u, :ubar) == 1.0
        @test get_flavor_factor(:sbar, :s) == 2.0
    end
    
    # 测试3：散射过程解析
    @testset "散射过程解析" begin
        # 标准格式
        @test parse_scattering_process(:uu_to_uu) == (:u, :u, :u, :u)
        @test parse_scattering_process(:us_to_us) == (:u, :s, :u, :s)
        @test parse_scattering_process(:ud_to_ud) == (:u, :d, :u, :d)
        @test parse_scattering_process(:dd_to_dd) == (:d, :d, :d, :d)
        @test parse_scattering_process(:ss_to_ss) == (:s, :s, :s, :s)
        @test parse_scattering_process(:ds_to_ds) == (:d, :s, :d, :s)
        
        # 反夸克格式
        @test parse_scattering_process(:uubar_to_uubar) == (:u, :ubar, :u, :ubar)
        @test parse_scattering_process(:ddbar_to_ddbar) == (:d, :dbar, :d, :dbar)
    end
    
    # 测试4：散射道味因子确定
    @testset "散射道味因子确定" begin
        # uu → uu 过程
        T1_t, T2_t = get_flavor_factors_for_channel(:uu_to_uu, :t)
        @test T1_t == 1.0  # u→u
        @test T2_t == 1.0  # u→u
        
        T1_u, T2_u = get_flavor_factors_for_channel(:uu_to_uu, :u)
        @test T1_u == 1.0  # u→u
        @test T2_u == 1.0  # u→u
        
        T1_s, T2_s = get_flavor_factors_for_channel(:uu_to_uu, :s)
        @test T1_s == 1.0  # u-u
        @test T2_s == 1.0  # u-u
        
        # us → us 过程
        T1_t, T2_t = get_flavor_factors_for_channel(:us_to_us, :t)
        @test T1_t == 1.0  # u→u
        @test T2_t == 2.0  # s→s
        
        T1_u, T2_u = get_flavor_factors_for_channel(:us_to_us, :u)
        @test T1_u ≈ √2   # u→s
        @test T2_u ≈ √2   # s→u
        
        # ud → ud 过程
        T1_u, T2_u = get_flavor_factors_for_channel(:ud_to_ud, :u)
        @test T1_u ≈ √2   # u→d
        @test T2_u ≈ √2   # d→u
    end
    
    # 测试5：质量提取功能（新增）
    @testset "质量提取功能" begin
        quark_params = (m = (u=0.3, d=0.3, s=0.5),)
        
        # qq散射过程
        m1, m2, m3, m4 = get_quark_masses_for_process(:uu_to_uu, quark_params)
        @test (m1, m2, m3, m4) == (0.3, 0.3, 0.3, 0.3)
        
        m1, m2, m3, m4 = get_quark_masses_for_process(:us_to_us, quark_params)
        @test (m1, m2, m3, m4) == (0.3, 0.5, 0.3, 0.5)
        
        m1, m2, m3, m4 = get_quark_masses_for_process(:ss_to_ss, quark_params)
        @test (m1, m2, m3, m4) == (0.5, 0.5, 0.5, 0.5)
        
        # qqbar散射过程（反夸克使用相同质量）
        m1, m2, m3, m4 = get_quark_masses_for_process(:uubar_to_uubar, quark_params)
        @test (m1, m2, m3, m4) == (0.3, 0.3, 0.3, 0.3)
        
        m1, m2, m3, m4 = get_quark_masses_for_process(:uubar_to_ssbar, quark_params)
        @test (m1, m2, m3, m4) == (0.3, 0.3, 0.5, 0.5)
    end
    
    # 测试6：质心系动量计算（新增）
    @testset "质心系动量计算" begin
        quark_params = (m = (u=0.3, d=0.3, s=0.5),)
        
        s = 4.0  # fm⁻²
        t = -0.5  # fm⁻²
        u = 4 * 0.3^2 - s - t  # 自动计算u
        
        # 测试s道（k=0）
        result = calculate_cms_momentum(:uu_to_uu, s, t, :s, quark_params)
        @test result.k == 0.0
        @test result.k0 > 0.0
        
        # 测试t道
        result_t = calculate_cms_momentum(:uu_to_uu, s, t, :t, quark_params)
        @test result_t.k0 > 0.0
        @test result_t.k >= 0.0
        
        # 测试u道
        result_u = calculate_cms_momentum(:uu_to_uu, s, t, :u, quark_params)
        @test result_u.k0 > 0.0
        @test result_u.k >= 0.0
        
        # 测试手动提供u参数
        result_manual_u = calculate_cms_momentum(:uu_to_uu, s, t, :u, quark_params; u=u)
        @test result_manual_u.k0 ≈ result_u.k0
        @test result_manual_u.k ≈ result_u.k
        
        # 测试边界条件（k0²-t略小于0的情况）
        # 设置一个使k0²-t接近0的参数
        s_boundary = 1.0
        t_boundary = -0.999  # 使得k0²-t可能为负
        result_boundary = calculate_cms_momentum(:uu_to_uu, s_boundary, t_boundary, :t, quark_params)
        @test result_boundary.k >= 0.0  # 应该设为0而不是NaN
    end
    
    # 测试7：通道分离接口（新增）
    @testset "通道分离接口" begin
        # 准备完整的物理参数
        T = 150.0 / 197.327  # 150 MeV
        m_u = 300.0 / 197.327
        m_s = 500.0 / 197.327
        μ_u = 0.0
        μ_s = 0.0
        Φ = 0.5
        Φbar = 0.5
        ξ = 0.0
        
        A_u = A(T, μ_u, m_u, Φ, Φbar)
        A_s = A(T, μ_s, m_s, Φ, Φbar)
        
        G_u = calculate_G_from_A(A_u)
        G_s = calculate_G_from_A(A_s)
        
        quark_params = (
            m = (u=m_u, d=m_u, s=m_s),
            μ = (u=μ_u, d=μ_u, s=μ_s),
            A = (u=A_u, d=A_u, s=A_s)
        )
        
        thermo_params = (T=T, Φ=Φ, Φbar=Φbar, ξ=ξ)
        
        # 使用真实的PNJL参数
        G_fm2 = Constants_PNJL.G_GeV_inv2 / (197.327^2)
        K_fm5 = Constants_PNJL.K_GeV_inv5 / (197.327^5)
        K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
        
        k0 = 100.0 / 197.327
        k_norm = 50.0 / 197.327
        
        # 测试qq散射（返回t_S, t_P, u_S, u_P）
        result_qq = calculate_all_propagators_by_channel(
            :uu_to_uu, k0, k_norm, quark_params, thermo_params, K_coeffs
        )
        @test haskey(result_qq, :t_S) && haskey(result_qq, :t_P)
        @test haskey(result_qq, :u_S) && haskey(result_qq, :u_P)
        @test result_qq.t_S isa ComplexF64
        @test result_qq.t_P isa ComplexF64
        @test result_qq.u_S isa ComplexF64
        @test result_qq.u_P isa ComplexF64
        
        # 测试qqbar散射（返回t_S, t_P, s_S, s_P）
        result_qqbar = calculate_all_propagators_by_channel(
            :uubar_to_uubar, k0, k_norm, quark_params, thermo_params, K_coeffs
        )
        @test haskey(result_qqbar, :t_S) && haskey(result_qqbar, :t_P)
        @test haskey(result_qqbar, :s_S) && haskey(result_qqbar, :s_P)
        @test result_qqbar.t_S isa ComplexF64
        @test result_qqbar.t_P isa ComplexF64
        @test result_qqbar.s_S isa ComplexF64
        @test result_qqbar.s_P isa ComplexF64
        
        # 验证S和P通道分离的一致性
        # D_total = D_S + D_P 应该与旧接口一致
        result_old = calculate_all_propagators(
            :uu_to_uu, k0, k_norm, quark_params, thermo_params, K_coeffs
        )
        D_t_total = result_qq.t_S + result_qq.t_P
        D_u_total = result_qq.u_S + result_qq.u_P
        @test D_t_total ≈ result_old.t
        @test D_u_total ≈ result_old.u
    end
    
    # 测试8：错误处理
    @testset "错误处理" begin
        quark_params = (m = (u=0.3, d=0.3, s=0.5),)
        
        # 无效的散射道
        @test_throws ErrorException calculate_cms_momentum(
            :uu_to_uu, 4.0, -0.5, :invalid, quark_params
        )
        
        # 无效的散射过程格式
        @test_throws ErrorException parse_scattering_process(:invalid_format)
        
        # 未知的散射过程
        T = 150.0 / 197.327
        k0 = 100.0 / 197.327
        k_norm = 50.0 / 197.327
        quark_params_full = (
            m = (u=0.3, d=0.3, s=0.5),
            μ = (u=0.0, d=0.0, s=0.0),
            A = (u=1.0, d=1.0, s=1.0)
        )
        thermo_params = (T=T, Φ=0.5, Φbar=0.5, ξ=0.0)
        K_coeffs = (K0_plus = 1.0, K0_minus = 0.8,
                    K123_plus = 0.9, K123_minus = 0.7,
                    K4567_plus = 0.85, K4567_minus = 0.65,
                    K8_plus = 0.95, K8_minus = 0.75,
                    K08_plus = 0.1, K08_minus = -0.1,
                    det_K_plus = 1.0, det_K_minus = 0.8)
        
        @test_throws ErrorException calculate_all_propagators(
            :unknown_process, k0, k_norm, quark_params_full, thermo_params, K_coeffs
        )
    end
    
    # 测试9：所有11种散射过程（新增）
    @testset "所有11种散射过程" begin
        T = 150.0 / 197.327
        m_u = 300.0 / 197.327
        m_s = 500.0 / 197.327
        μ_u = 0.0
        μ_s = 0.0
        Φ = 0.5
        Φbar = 0.5
        ξ = 0.0
        
        A_u = A(T, μ_u, m_u, Φ, Φbar)
        A_s = A(T, μ_s, m_s, Φ, Φbar)
        G_u = calculate_G_from_A(A_u)
        G_s = calculate_G_from_A(A_s)
        
        quark_params = (
            m = (u=m_u, d=m_u, s=m_s),
            μ = (u=μ_u, d=μ_u, s=μ_s),
            A = (u=A_u, d=A_u, s=A_s)
        )
        thermo_params = (T=T, Φ=Φ, Φbar=Φbar, ξ=ξ)
        
        G_fm2 = Constants_PNJL.G_GeV_inv2 / (197.327^2)
        K_fm5 = Constants_PNJL.K_GeV_inv5 / (197.327^5)
        K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
        
        k0 = 100.0 / 197.327
        k_norm = 50.0 / 197.327
        
        # 4种qq散射
        qq_processes = [:uu_to_uu, :ss_to_ss, :ud_to_ud, :us_to_us]
        for process in qq_processes
            result = calculate_all_propagators_by_channel(
                process, k0, k_norm, quark_params, thermo_params, K_coeffs
            )
            @test haskey(result, :t_S) && haskey(result, :t_P)
            @test haskey(result, :u_S) && haskey(result, :u_P)
            @test !isnan(real(result.t_S)) && !isnan(imag(result.t_S))
            @test !isnan(real(result.t_P)) && !isnan(imag(result.t_P))
            @test !isnan(real(result.u_S)) && !isnan(imag(result.u_S))
            @test !isnan(real(result.u_P)) && !isnan(imag(result.u_P))
        end
        
        # 7种qqbar散射
        qqbar_processes = [:udbar_to_udbar, :usbar_to_usbar, :uubar_to_uubar,
                          :uubar_to_ddbar, :uubar_to_ssbar, :ssbar_to_uubar, :ssbar_to_ssbar]
        for process in qqbar_processes
            result = calculate_all_propagators_by_channel(
                process, k0, k_norm, quark_params, thermo_params, K_coeffs
            )
            @test haskey(result, :t_S) && haskey(result, :t_P)
            @test haskey(result, :s_S) && haskey(result, :s_P)
            @test !isnan(real(result.t_S)) && !isnan(imag(result.t_S))
            @test !isnan(real(result.t_P)) && !isnan(imag(result.t_P))
            @test !isnan(real(result.s_S)) && !isnan(imag(result.s_S))
            @test !isnan(real(result.s_P)) && !isnan(imag(result.s_P))
        end
    end
end

println("\n" * "="^70)
println("总传播子计算测试完成！")
println("="^70)
