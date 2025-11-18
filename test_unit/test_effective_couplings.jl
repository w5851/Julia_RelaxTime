"""
测试 EffectiveCouplings 模块的有效耦合系数计算

运行方式：
```julia
using Pkg
Pkg.activate(".")
include("test/test_effective_couplings.jl")
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

using .Constants_PNJL: G_fm2, K_fm5, ħc_MeV_fm, N_color
using .GaussLegendre: gauleg
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings, coupling_matrix_determinant

# 标准测试参数
const T_TEST_MeV = 150.0
const μ_TEST_MeV = 0.0
const m_u_TEST_MeV = 300.0  # u/d 夸克有效质量
const m_s_TEST_MeV = 500.0  # s 夸克有效质量
const Φ_TEST = 0.5
const Φbar_TEST = 0.5

@testset "EffectiveCouplings 模块测试" verbose=true begin
    
    # 预生成积分节点（所有测试共用）
    nodes_p, weights_p = gauleg(0.0, 20.0, 64)
    
    @testset "基本功能测试" begin
        @testset "calculate_G_from_A 基本功能" begin
            # 测试典型A值
            A_typical = 1.5  # fm
            G_result = calculate_G_from_A(A_typical)
            
            @test G_result isa Float64
            @test isfinite(G_result)
            @test G_result < 0  # 应该是负值
            
            # 手动验证公式
            expected_G = -N_color / (4.0 * π^2) * A_typical
            @test abs(G_result - expected_G) < 1e-10
            
            @info "calculate_G_from_A 基本测试" A=A_typical G=G_result expected=expected_G
        end
        
        @testset "calculate_effective_couplings 基本功能" begin
            # 使用模拟的G^f值
            G_u_test = -0.3
            G_s_test = -0.2
            
            K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u_test, G_s_test)
            
            # 检查返回类型
            @test K_coeffs isa NamedTuple
            @test hasfield(typeof(K_coeffs), :K0_plus)
            @test hasfield(typeof(K_coeffs), :K123_minus)
            @test hasfield(typeof(K_coeffs), :K08_plus)
            
            # 检查所有系数都是有限数值
            for (key, val) in pairs(K_coeffs)
                @test isfinite(val)
                @test val isa Float64
            end
            
            @info "calculate_effective_couplings 基本测试" K_coeffs
        end
        
        @testset "coupling_matrix_determinant 基本功能" begin
            # 使用模拟值
            K0_plus_test = 1e-5
            K8_plus_test = 1.1e-5
            K08_plus_test = 1e-7
            
            det_K = coupling_matrix_determinant(K0_plus_test, K8_plus_test, K08_plus_test)
            
            @test det_K isa Float64
            @test isfinite(det_K)
            
            # 手动验证公式
            expected_det = K0_plus_test * K8_plus_test - K08_plus_test^2
            @test abs(det_K - expected_det) < 1e-20
            
            @info "coupling_matrix_determinant 基本测试" det_K=det_K expected=expected_det
        end
    end
    
    @testset "手征极限测试（K=0或G^f=0）" begin
        @testset "K=0 时所有K_α退化为G" begin
            G_u_test = -0.3
            G_s_test = -0.2
            K_zero = 0.0
            
            K_coeffs = calculate_effective_couplings(G_fm2, K_zero, G_u_test, G_s_test)
            
            # 当K=0时，所有K_α^±都应该等于G
            @test abs(K_coeffs.K0_plus - G_fm2) < 1e-15
            @test abs(K_coeffs.K0_minus - G_fm2) < 1e-15
            @test abs(K_coeffs.K123_plus - G_fm2) < 1e-15
            @test abs(K_coeffs.K123_minus - G_fm2) < 1e-15
            @test abs(K_coeffs.K4567_plus - G_fm2) < 1e-15
            @test abs(K_coeffs.K4567_minus - G_fm2) < 1e-15
            @test abs(K_coeffs.K8_plus - G_fm2) < 1e-15
            @test abs(K_coeffs.K8_minus - G_fm2) < 1e-15
            @test abs(K_coeffs.K08_plus) < 1e-15
            @test abs(K_coeffs.K08_minus) < 1e-15
            
            @info "K=0 手征极限测试通过" K_coeffs.K123_minus G_fm2
        end
        
        @testset "G^f=0 时所有K_α退化为G" begin
            G_u_zero = 0.0
            G_s_zero = 0.0
            
            K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u_zero, G_s_zero)
            
            # 当G^f=0时，所有K_α^±都应该等于G
            @test abs(K_coeffs.K0_plus - G_fm2) < 1e-15
            @test abs(K_coeffs.K0_minus - G_fm2) < 1e-15
            @test abs(K_coeffs.K123_plus - G_fm2) < 1e-15
            @test abs(K_coeffs.K123_minus - G_fm2) < 1e-15
            @test abs(K_coeffs.K8_plus - G_fm2) < 1e-15
            @test abs(K_coeffs.K8_minus - G_fm2) < 1e-15
            @test abs(K_coeffs.K08_plus) < 1e-15
            @test abs(K_coeffs.K08_minus) < 1e-15
            
            @info "G^f=0 手征极限测试通过"
        end
    end
    
    @testset "SU(3)对称性测试" begin
        @testset "G^u = G^s 时味道简并" begin
            G_symmetric = -0.25
            
            K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_symmetric, G_symmetric)
            
            # 当G^u = G^s时，π、K、η_8应该简并
            @test abs(K_coeffs.K123_plus - K_coeffs.K4567_plus) < 1e-15
            @test abs(K_coeffs.K123_minus - K_coeffs.K4567_minus) < 1e-15
            @test abs(K_coeffs.K123_plus - K_coeffs.K8_plus) < 1e-15
            @test abs(K_coeffs.K123_minus - K_coeffs.K8_minus) < 1e-15
            
            # 混合项应该为零
            @test abs(K_coeffs.K08_plus) < 1e-15
            @test abs(K_coeffs.K08_minus) < 1e-15
            
            @info "SU(3)对称性测试通过" G_u=G_symmetric K123=K_coeffs.K123_minus K4567=K_coeffs.K4567_minus
        end
    end
    
    @testset "物理约束测试" begin
        @testset "det K 正定性（标准参数）" begin
            # 使用真实物理参数计算
            T_inv_fm = T_TEST_MeV / ħc_MeV_fm
            μ_inv_fm = μ_TEST_MeV / ħc_MeV_fm
            m_u_inv_fm = m_u_TEST_MeV / ħc_MeV_fm
            m_s_inv_fm = m_s_TEST_MeV / ħc_MeV_fm
            
            # 计算A函数
            A_u = A(m_u_inv_fm, μ_inv_fm, T_inv_fm, Φ_TEST, Φbar_TEST, 
                    nodes_p, weights_p)
            A_s = A(m_s_inv_fm, μ_inv_fm, T_inv_fm, Φ_TEST, Φbar_TEST, 
                    nodes_p, weights_p)
            
            # 计算G^f
            G_u = calculate_G_from_A(A_u)
            G_s = calculate_G_from_A(A_s)
            
            # 计算K系数
            K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
            
            # 计算标量和赝标量通道的行列式
            det_K_scalar = coupling_matrix_determinant(K_coeffs.K0_plus, 
                                                       K_coeffs.K8_plus, 
                                                       K_coeffs.K08_plus)
            det_K_pseudoscalar = coupling_matrix_determinant(K_coeffs.K0_minus, 
                                                             K_coeffs.K8_minus, 
                                                             K_coeffs.K08_minus)
            
            @test det_K_scalar > 0
            @test det_K_pseudoscalar > 0
            
            @info "det K 正定性测试" scalar=det_K_scalar pseudoscalar=det_K_pseudoscalar
        end
        
        @testset "K系数量级检验" begin
            G_u_test = -0.3
            G_s_test = -0.2
            
            K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u_test, G_s_test)
            
            # K系数应该在合理范围内（与G同量级）
            for (key, val) in pairs(K_coeffs)
                if key != :K08_plus && key != :K08_minus
                    # 主要通道的K系数应该接近G（相对偏差<50%）
                    rel_diff = abs(val - G_fm2) / G_fm2
                    @test rel_diff < 0.5
                else
                    # 混合项应该远小于主要项
                    @test abs(val) < 0.1 * G_fm2
                end
            end
            
            @info "K系数量级检验通过"
        end
    end
    
    @testset "数值稳定性测试" begin
        @testset "极端温度扫描" begin
            T_values_MeV = [50.0, 100.0, 150.0, 200.0, 250.0]
            
            println("\n" * "="^80)
            println("温度扫描测试（μ=0）")
            println("="^80)
            println(@sprintf("%-12s %-12s %-12s %-12s %-12s", 
                           "T (MeV)", "G^u", "K₁₂₃⁻ (fm²)", "K₁₂₃⁻/G", "det K^P (fm⁴)"))
            println("-"^80)
            
            for T_MeV in T_values_MeV
                T_inv_fm = T_MeV / ħc_MeV_fm
                μ_inv_fm = 0.0
                m_u_inv_fm = m_u_TEST_MeV / ħc_MeV_fm
                m_s_inv_fm = m_s_TEST_MeV / ħc_MeV_fm
                
                # 计算A和G^f
                A_u = A(m_u_inv_fm, μ_inv_fm, T_inv_fm, 0.5, 0.5, nodes_p, weights_p)
                A_s = A(m_s_inv_fm, μ_inv_fm, T_inv_fm, 0.5, 0.5, nodes_p, weights_p)
                G_u = calculate_G_from_A(A_u)
                G_s = calculate_G_from_A(A_s)
                
                # 计算K系数
                K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
                
                # 计算行列式
                det_K = coupling_matrix_determinant(K_coeffs.K0_minus, 
                                                    K_coeffs.K8_minus, 
                                                    K_coeffs.K08_minus)
                
                # 检查数值合理性
                @test isfinite(G_u)
                @test isfinite(K_coeffs.K123_minus)
                @test det_K > 0
                
                ratio = K_coeffs.K123_minus / G_fm2
                println(@sprintf("%-12.1f %-12.4f %-12.3e %-12.4f %-12.3e", 
                               T_MeV, G_u, K_coeffs.K123_minus, ratio, det_K))
            end
            println()
        end
        
        @testset "化学势扫描" begin
            μ_values_MeV = [0.0, 100.0, 200.0, 300.0]
            
            println("\n" * "="^80)
            println("化学势扫描测试（T=150 MeV）")
            println("="^80)
            println(@sprintf("%-12s %-12s %-12s %-12s %-12s", 
                           "μ (MeV)", "G^u", "G^s", "|K₀₈⁻| (fm²)", "det K^P (fm⁴)"))
            println("-"^80)
            
            T_inv_fm = 150.0 / ħc_MeV_fm
            
            for μ_MeV in μ_values_MeV
                μ_inv_fm = μ_MeV / ħc_MeV_fm
                m_u_inv_fm = m_u_TEST_MeV / ħc_MeV_fm
                m_s_inv_fm = m_s_TEST_MeV / ħc_MeV_fm
                
                # 计算A和G^f
                A_u = A(m_u_inv_fm, μ_inv_fm, T_inv_fm, 0.5, 0.5, nodes_p, weights_p)
                A_s = A(m_s_inv_fm, μ_inv_fm, T_inv_fm, 0.5, 0.5, nodes_p, weights_p)
                G_u = calculate_G_from_A(A_u)
                G_s = calculate_G_from_A(A_s)
                
                # 计算K系数
                K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
                
                # 计算行列式
                det_K = coupling_matrix_determinant(K_coeffs.K0_minus, 
                                                    K_coeffs.K8_minus, 
                                                    K_coeffs.K08_minus)
                
                # 检查数值合理性
                @test isfinite(G_u)
                @test isfinite(G_s)
                @test det_K > 0
                
                println(@sprintf("%-12.1f %-12.4f %-12.4f %-12.3e %-12.3e", 
                               μ_MeV, G_u, G_s, abs(K_coeffs.K08_minus), det_K))
            end
            println()
        end
    end
    
    @testset "完整计算流程演示" begin
        println("\n" * "="^80)
        println("完整计算流程：从A函数到有效耦合系数")
        println("="^80)
        
        # 物理参数
        T_MeV = T_TEST_MeV
        μ_MeV = μ_TEST_MeV
        m_u_MeV = m_u_TEST_MeV
        m_s_MeV = m_s_TEST_MeV
        
        T_inv_fm = T_MeV / ħc_MeV_fm
        μ_inv_fm = μ_MeV / ħc_MeV_fm
        m_u_inv_fm = m_u_MeV / ħc_MeV_fm
        m_s_inv_fm = m_s_MeV / ħc_MeV_fm
        
        println("输入参数:")
        println("  T = $T_MeV MeV")
        println("  μ = $μ_MeV MeV")
        println("  m_u = $m_u_MeV MeV")
        println("  m_s = $m_s_MeV MeV")
        println("  Φ = $(Φ_TEST), Φ̄ = $(Φbar_TEST)")
        println()
        
        # 步骤1: 计算A函数
        println("步骤1: 计算单圈积分A")
        A_u = A(m_u_inv_fm, μ_inv_fm, T_inv_fm, Φ_TEST, Φbar_TEST, 
                nodes_p, weights_p)
        A_s = A(m_s_inv_fm, μ_inv_fm, T_inv_fm, Φ_TEST, Φbar_TEST, 
                nodes_p, weights_p)
        println("  A_u = $(@sprintf("%.6f", A_u)) fm")
        println("  A_s = $(@sprintf("%.6f", A_s)) fm")
        println()
        
        # 步骤2: 计算G^f
        println("步骤2: 计算夸克凝聚函数")
        G_u = calculate_G_from_A(A_u)
        G_s = calculate_G_from_A(A_s)
        println("  G^u = $(@sprintf("%.6f", G_u))")
        println("  G^s = $(@sprintf("%.6f", G_s))")
        println()
        
        # 步骤3: 计算K系数
        println("步骤3: 计算有效耦合系数")
        K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
        println("  K₀⁺   = $(@sprintf("%.3e", K_coeffs.K0_plus)) fm² (单态标量)")
        println("  K₀⁻   = $(@sprintf("%.3e", K_coeffs.K0_minus)) fm² (单态赝标量)")
        println("  K₁₂₃⁺ = $(@sprintf("%.3e", K_coeffs.K123_plus)) fm² (π通道标量)")
        println("  K₁₂₃⁻ = $(@sprintf("%.3e", K_coeffs.K123_minus)) fm² (π介子)")
        println("  K₄₅₆₇⁺ = $(@sprintf("%.3e", K_coeffs.K4567_plus)) fm² (K通道标量)")
        println("  K₄₅₆₇⁻ = $(@sprintf("%.3e", K_coeffs.K4567_minus)) fm² (K介子)")
        println("  K₈⁺   = $(@sprintf("%.3e", K_coeffs.K8_plus)) fm² (八重态标量)")
        println("  K₈⁻   = $(@sprintf("%.3e", K_coeffs.K8_minus)) fm² (八重态赝标量)")
        println("  K₀₈⁺  = $(@sprintf("%.3e", K_coeffs.K08_plus)) fm² (混合标量)")
        println("  K₀₈⁻  = $(@sprintf("%.3e", K_coeffs.K08_minus)) fm² (η-η'混合)")
        println()
        
        # 步骤4: 计算行列式
        println("步骤4: 检验物理约束")
        det_K_scalar = coupling_matrix_determinant(K_coeffs.K0_plus, 
                                                   K_coeffs.K8_plus, 
                                                   K_coeffs.K08_plus)
        det_K_pseudoscalar = coupling_matrix_determinant(K_coeffs.K0_minus, 
                                                         K_coeffs.K8_minus, 
                                                         K_coeffs.K08_minus)
        println("  det K^S = $(@sprintf("%.3e", det_K_scalar)) fm⁴")
        println("  det K^P = $(@sprintf("%.3e", det_K_pseudoscalar)) fm⁴")
        
        if det_K_scalar > 0 && det_K_pseudoscalar > 0
            println("  ✅ 物理约束满足：det K > 0")
        else
            println("  ⚠️ 警告：det K ≤ 0，模型可能失效")
        end
        println("="^80)
        println()
        
        @test det_K_scalar > 0
        @test det_K_pseudoscalar > 0
    end
end

println("\n所有测试完成！")
