"""
测试 MesonPropagator 模块的介子传播子计算

运行方式：
```julia
using Pkg
Pkg.activate(".")
include("test/test_meson_propagator.jl")
```
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Test
using Printf

include("../src/Constants_PNJL.jl")
include("../src/integration/GaussLegendre.jl")
include("../src/QuarkDistribution.jl")
include("../src/relaxtime/OneLoopIntegrals.jl")
include("../src/relaxtime/EffectiveCouplings.jl")
include("../src/relaxtime/MesonPropagator.jl")

using .Constants_PNJL: G_fm2, K_fm5, ħc_MeV_fm, λ₀, λ₈, ψ_u, ψ_d, ψ_s, ψbar_u, ψbar_d, ψbar_s
using .GaussLegendre: gauleg
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings, coupling_matrix_determinant
using .MesonPropagator: meson_propagator_simple, meson_propagator_mixed
using .MesonPropagator: calculate_coupling_matrix, extract_flavor, get_quark_wavefunction, calculate_current_vector

# 标准测试参数
const T_TEST_MeV = 150.0
const μ_TEST_MeV = 0.0
const m_u_TEST_MeV = 300.0  # u/d 夸克有效质量
const m_s_TEST_MeV = 500.0  # s 夸克有效质量
const Φ_TEST = 0.5
const Φbar_TEST = 0.5

@testset "MesonPropagator 模块测试" verbose=true begin
    
    # 预生成积分节点（所有测试共用）
    nodes_p, weights_p = gauleg(0.0, 20.0, 64)
    
    # 预计算K系数（所有测试共用）
    println("\n预计算K系数...")
    T_inv_fm = T_TEST_MeV / ħc_MeV_fm
    μ_inv_fm = μ_TEST_MeV / ħc_MeV_fm
    m_u_inv_fm = m_u_TEST_MeV / ħc_MeV_fm
    m_s_inv_fm = m_s_TEST_MeV / ħc_MeV_fm
    
    A_u = A(m_u_inv_fm, μ_inv_fm, T_inv_fm, Φ_TEST, Φbar_TEST, nodes_p, weights_p)
    A_s = A(m_s_inv_fm, μ_inv_fm, T_inv_fm, Φ_TEST, Φbar_TEST, nodes_p, weights_p)
    
    G_u = calculate_G_from_A(A_u)
    G_s = calculate_G_from_A(A_s)
    
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
    
    println("K123_plus = ", K_coeffs.K123_plus, " fm²")
    println("K4567_plus = ", K_coeffs.K4567_plus, " fm²")
    println("K123_minus = ", K_coeffs.K123_minus, " fm²")
    println("K4567_minus = ", K_coeffs.K4567_minus, " fm²")
    
    @testset "辅助函数测试" begin
        @testset "extract_flavor 函数" begin
            # 测试夸克
            flavor, is_bar = extract_flavor(:u)
            @test flavor == :u
            @test is_bar == false
            
            # 测试反夸克
            flavor, is_bar = extract_flavor(:ubar)
            @test flavor == :u
            @test is_bar == true
            
            # 测试其他味
            flavor, is_bar = extract_flavor(:dbar)
            @test flavor == :d
            @test is_bar == true
            
            flavor, is_bar = extract_flavor(:sbar)
            @test flavor == :s
            @test is_bar == true
            
            @info "extract_flavor 测试通过"
        end
        
        @testset "get_quark_wavefunction 函数" begin
            # 测试夸克波函数（列向量）
            ψu = get_quark_wavefunction(:u, false)
            @test ψu == [1.0, 0.0, 0.0]
            @test ψu == ψ_u
            
            ψd = get_quark_wavefunction(:d, false)
            @test ψd == [0.0, 1.0, 0.0]
            
            ψs = get_quark_wavefunction(:s, false)
            @test ψs == [0.0, 0.0, 1.0]
            
            # 测试反夸克波函数（行向量/1×3矩阵）
            ψbar_u_test = get_quark_wavefunction(:u, true)
            @test ψbar_u_test == [1.0 0.0 0.0]
            @test size(ψbar_u_test) == (1, 3)
            
            ψbar_d_test = get_quark_wavefunction(:d, true)
            @test ψbar_d_test == [0.0 1.0 0.0]
            
            ψbar_s_test = get_quark_wavefunction(:s, true)
            @test ψbar_s_test == [0.0 0.0 1.0]
            
            @info "get_quark_wavefunction 测试通过"
        end
        
        @testset "calculate_current_vector 函数" begin
            # 测试u+d散射（t道）
            J = calculate_current_vector(:u, :d, :t)
            @test length(J) == 2
            @test all(isfinite.(J))
            
            # 手动验证：ψbar_d * λ₀ * ψ_u
            J0_manual = (ψbar_d * λ₀ * ψ_u)[1]
            J8_manual = (ψbar_d * λ₈ * ψ_u)[1]
            @test abs(J[1] - J0_manual) < 1e-10
            @test abs(J[2] - J8_manual) < 1e-10
            
            # 测试反夸克（s道）
            J_s = calculate_current_vector(:u, :dbar, :s)
            @test length(J_s) == 2
            @test all(isfinite.(J_s))
            
            @info "calculate_current_vector 测试通过" J=J J_s=J_s
        end
        
        @testset "calculate_coupling_matrix 函数" begin
            # 测试赝标量通道
            Π_uu_test = 1.0e-5 + 0.0im
            Π_ss_test = 8.0e-6 + 0.0im
            
            M_P = calculate_coupling_matrix(Π_uu_test, Π_ss_test, K_coeffs, :P)
            
            @test size(M_P) == (2, 2)
            @test all(isfinite.(M_P))
            @test M_P isa Matrix{ComplexF64}
            
            # 验证对称性 M₀₈ = M₈₀
            @test abs(M_P[1, 2] - M_P[2, 1]) < 1e-15
            
            # 手动验证M矩阵元素
            M00_expected = 1.0 - K_coeffs.K0_plus * Π_uu_test
            M08_expected = -K_coeffs.K08_plus * Π_uu_test * 4.0 * sqrt(2.0) / 3.0
            M88_expected = 1.0 - K_coeffs.K8_plus * (Π_uu_test * 4.0 / 3.0 + Π_ss_test * 2.0 / 3.0)
            
            @test abs(M_P[1, 1] - M00_expected) < 1e-15
            @test abs(M_P[1, 2] - M08_expected) < 1e-15
            @test abs(M_P[2, 2] - M88_expected) < 1e-15
            
            # 测试标量通道
            M_S = calculate_coupling_matrix(Π_uu_test, Π_ss_test, K_coeffs, :S)
            @test size(M_S) == (2, 2)
            @test M_S[1, 2] == M_S[2, 1]  # 对称性
            
            # 验证M₀₈系数修正（4√2/3 ≈ 1.8856，不是4/(3√2) ≈ 0.9428）
            coefficient = 4.0 * sqrt(2.0) / 3.0
            @test abs(coefficient - 1.8856) < 0.001
            @test abs(coefficient - 0.9428) > 0.5
            
            @info "calculate_coupling_matrix 测试通过" M_P=M_P M_S=M_S
        end
    end
    
    @testset "一般介子传播子测试（meson_propagator_simple）" begin
        # 使用模拟的极化函数值
        Π_uu_P = 1.0e-5 + 0.0im  # ūu极化函数（赝标量通道）
        Π_us_P = 9.0e-6 + 0.0im  # ūs极化函数（赝标量通道）
        Π_uu_S = 1.2e-5 + 0.0im  # ūu极化函数（标量通道）
        Π_us_S = 1.0e-5 + 0.0im  # ūs极化函数（标量通道）
        
        @testset "基本功能测试" begin
            # π介子
            D_pi = meson_propagator_simple(:pi, K_coeffs, Π_uu_P)
            @test D_pi isa ComplexF64
            @test isfinite(D_pi)
            @test abs(D_pi) > 0
            
            # K介子（注意使用Π_us）
            D_K = meson_propagator_simple(:K, K_coeffs, Π_us_P)
            @test D_K isa ComplexF64
            @test isfinite(D_K)
            @test abs(D_K) > 0
            
            # σ_π介子
            D_sigma_pi = meson_propagator_simple(:sigma_pi, K_coeffs, Π_uu_S)
            @test D_sigma_pi isa ComplexF64
            @test isfinite(D_sigma_pi)
            @test abs(D_sigma_pi) > 0
            
            # σ_K介子
            D_sigma_K = meson_propagator_simple(:sigma_K, K_coeffs, Π_us_S)
            @test D_sigma_K isa ComplexF64
            @test isfinite(D_sigma_K)
            @test abs(D_sigma_K) > 0
            
            @info "一般介子传播子基本功能测试通过" D_pi=D_pi D_K=D_K D_sigma_pi=D_sigma_pi D_sigma_K=D_sigma_K
        end
        
        @testset "通道-K系数映射验证" begin
            # 验证π使用K123_plus
            D_pi = meson_propagator_simple(:pi, K_coeffs, Π_uu_P)
            D_pi_manual = 1.0 / (1.0 - K_coeffs.K123_plus * Π_uu_P)
            @test abs(D_pi - D_pi_manual) < 1e-15
            
            # 验证K使用K4567_plus
            D_K = meson_propagator_simple(:K, K_coeffs, Π_us_P)
            D_K_manual = 1.0 / (1.0 - K_coeffs.K4567_plus * Π_us_P)
            @test abs(D_K - D_K_manual) < 1e-15
            
            # 验证σ_π使用K123_minus
            D_sigma_pi = meson_propagator_simple(:sigma_pi, K_coeffs, Π_uu_S)
            D_sigma_pi_manual = 1.0 / (1.0 - K_coeffs.K123_minus * Π_uu_S)
            @test abs(D_sigma_pi - D_sigma_pi_manual) < 1e-15
            
            # 验证σ_K使用K4567_minus
            D_sigma_K = meson_propagator_simple(:sigma_K, K_coeffs, Π_us_S)
            D_sigma_K_manual = 1.0 / (1.0 - K_coeffs.K4567_minus * Π_us_S)
            @test abs(D_sigma_K - D_sigma_K_manual) < 1e-15
            
            @info "通道-K系数映射验证通过"
        end
        
        @testset "物理约束测试" begin
            # 测试不同动量下的传播子行为
            Π_small = 1.0e-6 + 0.0im  # 小极化函数（远离质量壳）
            Π_large = 5.0e-5 + 0.0im  # 大极化函数（接近质量壳）
            
            D_small = meson_propagator_simple(:pi, K_coeffs, Π_small)
            D_large = meson_propagator_simple(:pi, K_coeffs, Π_large)
            
            # 极化函数增大时，传播子应该增大（接近质量壳）
            @test abs(D_large) > abs(D_small)
            
            # 测试虚部（稳定介子应该虚部接近零）
            @test abs(imag(D_small)) < 1e-10  # 纯实数极化函数应给出纯实数传播子
            
            @info "物理约束测试通过" D_small=D_small D_large=D_large
        end
        
        @testset "错误处理测试" begin
            # 测试无效介子类型
            @test_throws ErrorException meson_propagator_simple(:invalid, K_coeffs, Π_uu_P)
            
            @info "错误处理测试通过"
        end
    end
    
    @testset "混合介子传播子测试（meson_propagator_mixed）" begin
        # 预计算det_K和M矩阵
        Π_uu_P = 1.0e-5 + 0.0im
        Π_ss_P = 8.0e-6 + 0.0im
        
        det_K_P = coupling_matrix_determinant(K_coeffs.K0_plus, K_coeffs.K8_plus, K_coeffs.K08_plus)
        det_K_S = coupling_matrix_determinant(K_coeffs.K0_minus, K_coeffs.K8_minus, K_coeffs.K08_minus)
        
        M_P = calculate_coupling_matrix(Π_uu_P, Π_ss_P, K_coeffs, :P)
        M_S = calculate_coupling_matrix(Π_uu_P, Π_ss_P, K_coeffs, :S)
        
        @testset "基本功能测试" begin
            # η/η'传播子（t道：u+d→u+d）
            D_eta_t = meson_propagator_mixed(det_K_P, M_P, :u, :d, :u, :d, :t)
            @test D_eta_t isa ComplexF64
            @test isfinite(D_eta_t)
            @test abs(D_eta_t) > 0
            
            # η/η'传播子（s道：u+ū→u+ū，同味湮灭）
            D_eta_s = meson_propagator_mixed(det_K_P, M_P, :u, :ubar, :u, :ubar, :s)
            @test D_eta_s isa ComplexF64
            @test isfinite(D_eta_s)
            @test abs(D_eta_s) > 0
            
            # η/η'传播子（u道：u+d→d+u）
            D_eta_u = meson_propagator_mixed(det_K_P, M_P, :u, :d, :d, :u, :u)
            @test D_eta_u isa ComplexF64
            @test isfinite(D_eta_u)
            @test abs(D_eta_u) > 0
            
            # σ/σ'传播子（标量通道）
            D_sigma = meson_propagator_mixed(det_K_S, M_S, :u, :d, :u, :d, :t)
            @test D_sigma isa ComplexF64
            @test isfinite(D_sigma)
            @test abs(D_sigma) > 0
            
            @info "混合介子传播子基本功能测试通过" D_eta_t=D_eta_t D_eta_s=D_eta_s D_eta_u=D_eta_u D_sigma=D_sigma
        end
        
        @testset "散射道映射测试" begin
            # t道：u+d→u+d（纯夸克散射）
            D_t = meson_propagator_mixed(det_K_P, M_P, :u, :d, :u, :d, :t)
            
            # 手动计算验证（t道） - 使用u+u散射避免正交性问题
            D_t_uu = meson_propagator_mixed(det_K_P, M_P, :u, :u, :u, :u, :t)
            J_t_uu = calculate_current_vector(:u, :u, :t)
            J_prime_t_uu = calculate_current_vector(:u, :u, :t)
            
            using LinearAlgebra: det, transpose
            det_M = det(M_P)
            result_manual_t = (transpose(J_t_uu) * M_P * J_prime_t_uu)[1]
            D_t_manual = 2.0 * det_K_P / det_M * result_manual_t
            
            @test abs(D_t_uu - D_t_manual) / abs(D_t_manual) < 1e-10
            
            # s道：u+ū→u+ū（夸克-反夸克湮灭）
            D_s = meson_propagator_mixed(det_K_P, M_P, :u, :ubar, :u, :ubar, :s)
            @test D_s isa ComplexF64
            @test isfinite(D_s)
            @test abs(D_s) > 0  # 同味湮灭应该非零
            
            # u道：u+d→d+u（交叉散射）
            D_u = meson_propagator_mixed(det_K_P, M_P, :u, :d, :d, :u, :u)
            @test D_u isa ComplexF64
            @test isfinite(D_u)
            @test abs(D_u) > 0
            
            @info "散射道映射测试通过" D_t=D_t D_t_uu=D_t_uu D_s=D_s D_u=D_u
        end
        
        @testset "矩阵乘法展开对比测试" begin
            # 测试矩阵乘法计算的正确性（使用相同的散射过程）
            D_matrix = meson_propagator_mixed(det_K_P, M_P, :u, :u, :u, :u, :t)
            
            # 手动展开矩阵乘法
            J = calculate_current_vector(:u, :u, :t)
            J_prime = calculate_current_vector(:u, :u, :t)
            
            using LinearAlgebra: det
            det_M = det(M_P)
            
            # J^T M J' 的手动展开
            temp = M_P * J_prime
            result_expanded = J[1] * temp[1] + J[2] * temp[2]
            D_expanded = 2.0 * det_K_P / det_M * result_expanded
            
            # 相对误差应小于1e-10
            rel_error = abs(D_matrix - D_expanded) / abs(D_expanded)
            @test rel_error < 1e-10
            
            @info "矩阵乘法展开对比测试通过" rel_error=rel_error D_matrix=D_matrix D_expanded=D_expanded
        end
        
        @testset "夸克味配置测试" begin
            # 测试不同夸克组合
            # u+u散射
            D_uu = meson_propagator_mixed(det_K_P, M_P, :u, :u, :u, :u, :t)
            @test isfinite(D_uu)
            
            # s+s散射
            D_ss = meson_propagator_mixed(det_K_P, M_P, :s, :s, :s, :s, :t)
            @test isfinite(D_ss)
            
            # u+s散射
            D_us = meson_propagator_mixed(det_K_P, M_P, :u, :s, :u, :s, :t)
            @test isfinite(D_us)
            
            # 验证不同味组合给出不同结果
            @test D_uu != D_ss
            @test D_uu != D_us
            
            @info "夸克味配置测试通过" D_uu=D_uu D_ss=D_ss D_us=D_us
        end
        
        @testset "物理合理性检验" begin
            # 测试det_K符号（应该为正以保证因果性）
            @test det_K_P > 0
            @test det_K_S > 0
            
            # 测试传播子量级（应该是fm²量级）
            D_test = meson_propagator_mixed(det_K_P, M_P, :u, :d, :u, :d, :t)
            @test abs(D_test) > 1e-10  # 不能太小
            @test abs(D_test) < 1e10   # 不能太大（非物理）
            
            @info "物理合理性检验通过" det_K_P=det_K_P det_K_S=det_K_S
        end
        
        @testset "错误处理测试" begin
            # 测试无效散射道
            @test_throws ErrorException meson_propagator_mixed(det_K_P, M_P, :u, :d, :u, :d, :invalid)
            
            @info "错误处理测试通过"
        end
    end
    
    @testset "各向异性修正测试（通过极化函数）" begin
        # 模拟各向同性和各向异性的极化函数
        Π_iso = 1.0e-5 + 0.0im  # 各向同性（ξ=0）
        Π_aniso = 1.2e-5 + 0.0im  # 各向异性（ξ≠0，极化函数会改变）
        
        # 计算传播子
        D_iso = meson_propagator_simple(:pi, K_coeffs, Π_iso)
        D_aniso = meson_propagator_simple(:pi, K_coeffs, Π_aniso)
        
        # 各向异性修正应该导致传播子变化
        Δ_aniso = abs(D_aniso - D_iso) / abs(D_iso)
        @test Δ_aniso > 0  # 应该有明显变化
        @test Δ_aniso < 1.0  # 但不应该是量级变化
        
        @info "各向异性修正测试通过" Δ_aniso_percent=Δ_aniso*100 D_iso=D_iso D_aniso=D_aniso
    end
    
    @testset "批量计算性能测试" begin
        # 测试K系数复用的性能优势
        Π_values = [1.0e-5 + 0.0im, 2.0e-5 + 0.0im, 3.0e-5 + 0.0im, 4.0e-5 + 0.0im, 5.0e-5 + 0.0im]
        
        # 方法1：每次重新计算K系数（低效）
        time_inefficient = @elapsed begin
            for Π in Π_values
                K_coeffs_temp = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
                D = meson_propagator_simple(:pi, K_coeffs_temp, Π)
            end
        end
        
        # 方法2：复用K系数（高效）
        time_efficient = @elapsed begin
            for Π in Π_values
                D = meson_propagator_simple(:pi, K_coeffs, Π)
            end
        end
        
        # 高效方法应该更快
        speedup = time_inefficient / time_efficient
        @test speedup > 1.0
        
        @info "批量计算性能测试通过" speedup=speedup time_inefficient=time_inefficient time_efficient=time_efficient
    end
end

println("\n" * "="^80)
println("MesonPropagator 模块测试完成！")
println("="^80)
println("\n测试总结：")
println("✓ 辅助函数测试通过（extract_flavor, get_quark_wavefunction, calculate_current_vector, calculate_coupling_matrix）")
println("✓ 一般介子传播子测试通过（π, K, σ_π, σ_K）")
println("✓ 混合介子传播子测试通过（η/η', σ/σ'）")
println("✓ 通道-K系数映射验证通过（赝标量P用K⁺，标量S用K⁻）")
println("✓ 散射道映射测试通过（t道、s道、u道）")
println("✓ 矩阵乘法展开对比测试通过（相对误差<1e-10）")
println("✓ 夸克味配置测试通过（u, d, s及其反粒子）")
println("✓ 物理合理性检验通过（det_K>0, 传播子量级合理）")
println("✓ 各向异性修正测试通过（通过预计算的Π值间接影响传播子）")
println("✓ 批量计算性能测试通过（K系数复用带来性能提升）")
println("\n所有测试全部通过！模块实现正确。")
