"""
测试总传播子计算模块

测试内容：
1. 味因子查询功能
2. 散射过程解析
3. 味因子自动确定（t/u/s道）
4. 一般介子总传播子计算
5. 混合介子总传播子计算
6. 复杂散射过程组合
"""

using Test

push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../src/relaxtime"))

include("../src/relaxtime/TotalPropagator.jl")
using .TotalPropagator
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
    
    # 测试5：一般介子总传播子计算
    @testset "一般介子总传播子" begin
        # 准备测试数据
        include("../src/relaxtime/EffectiveCouplings.jl")
        using .EffectiveCouplings
        
        # 模拟K系数
        K_coeffs = (
            K0_plus = 1.0, K0_minus = 0.8,
            K123_plus = 0.9, K123_minus = 0.7,
            K4567_plus = 0.85, K4567_minus = 0.65,
            K8_plus = 0.95, K8_minus = 0.75,
            K08_plus = 0.1, K08_minus = -0.1
        )
        
        # 模拟极化函数
        Π_pi = ComplexF64(0.1, 0.01)
        Π_K = ComplexF64(0.08, 0.008)
        Π_dict = Dict(:pi => Π_pi, :K => Π_K)
        
        # 计算uu→uu过程的t道传播子（π+K交换）
        D_total = total_propagator_simple(:uu_to_uu, :t, [:pi, :K], K_coeffs, Π_dict)
        
        # 验证结果是复数
        @test D_total isa ComplexF64
        
        # 手动计算验证
        include("../src/relaxtime/MesonPropagator.jl")
        using .MesonPropagator
        
        D_pi = meson_propagator_simple(:pi, K_coeffs, Π_pi)
        D_K = meson_propagator_simple(:K, K_coeffs, Π_K)
        T1, T2 = get_flavor_factors_for_channel(:uu_to_uu, :t)
        D_expected = T1 * (D_pi + D_K) * T2
        
        @test D_total ≈ D_expected
    end
    
    # 测试6：味因子的正确应用
    @testset "味因子应用验证" begin
        K_coeffs = (
            K0_plus = 1.0, K0_minus = 0.8,
            K123_plus = 0.9, K123_minus = 0.7,
            K4567_plus = 0.85, K4567_minus = 0.65,
            K8_plus = 0.95, K8_minus = 0.75,
            K08_plus = 0.1, K08_minus = -0.1
        )
        
        Π_pi = ComplexF64(0.1, 0.01)
        Π_dict = Dict(:pi => Π_pi)
        
        # us→us过程的t道：T1=1, T2=2
        D_us_t = total_propagator_simple(:us_to_us, :t, [:pi], K_coeffs, Π_dict)
        
        # us→us过程的u道：T1=√2, T2=√2
        D_us_u = total_propagator_simple(:us_to_us, :u, [:pi], K_coeffs, Π_dict)
        
        # 计算基础传播子
        include("../src/relaxtime/MesonPropagator.jl")
        using .MesonPropagator
        D_pi = meson_propagator_simple(:pi, K_coeffs, Π_pi)
        
        # 验证味因子的影响
        @test D_us_t ≈ 1.0 * D_pi * 2.0  # t道: T1=1, T2=2
        @test D_us_u ≈ √2 * D_pi * √2    # u道: T1=√2, T2=√2
        @test real(D_us_u) ≈ 2.0 * real(D_pi)  # √2 * √2 = 2
    end
    
    # 测试7：错误处理
    @testset "错误处理" begin
        K_coeffs = (K0_plus = 1.0, K0_minus = 0.8,
                    K123_plus = 0.9, K123_minus = 0.7,
                    K4567_plus = 0.85, K4567_minus = 0.65,
                    K8_plus = 0.95, K8_minus = 0.75,
                    K08_plus = 0.1, K08_minus = -0.1)
        
        # 极化函数字典缺少介子
        Π_dict_incomplete = Dict(:pi => ComplexF64(0.1, 0.01))
        
        @test_throws ErrorException total_propagator_simple(
            :uu_to_uu, :t, [:pi, :K], K_coeffs, Π_dict_incomplete
        )
        
        # 无效的散射道
        Π_dict_complete = Dict(:pi => ComplexF64(0.1, 0.01))
        @test_throws ErrorException total_propagator_simple(
            :uu_to_uu, :invalid, [:pi], K_coeffs, Π_dict_complete
        )
        
        # 无效的散射过程格式
        @test_throws ErrorException parse_scattering_process(:invalid_format)
    end
    
    # 测试8：复杂散射过程
    @testset "复杂散射过程" begin
        # 测试所有9种夸克组合
        processes = [:uu_to_uu, :ud_to_ud, :us_to_us,
                     :dd_to_dd, :ds_to_ds, :ss_to_ss,
                     :du_to_du, :su_to_su, :sd_to_sd]
        
        K_coeffs = (K0_plus = 1.0, K0_minus = 0.8,
                    K123_plus = 0.9, K123_minus = 0.7,
                    K4567_plus = 0.85, K4567_minus = 0.65,
                    K8_plus = 0.95, K8_minus = 0.75,
                    K08_plus = 0.1, K08_minus = -0.1)
        
        Π_dict = Dict(:pi => ComplexF64(0.1, 0.01))
        
        # 确保所有过程都能正确计算
        for process in processes
            for channel in [:t, :u, :s]
                D = total_propagator_simple(process, channel, [:pi], K_coeffs, Π_dict)
                @test D isa ComplexF64
                @test !isnan(real(D)) && !isnan(imag(D))
            end
        end
    end
end

println("\n" * "="^70)
println("总传播子计算测试完成！")
println("="^70)
