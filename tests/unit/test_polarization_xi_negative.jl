"""
诊断ξ<0时极化函数问题

测试内容：
1. 检查correction_cos_theta_coefficient对ξ符号的响应
2. 检查B0_correction对ξ符号的响应
3. 检查polarization_aniso对ξ符号的响应
4. 定位问题根源
"""

using Test
using Printf

push!(LOAD_PATH, joinpath(@__DIR__, "../../src"))

include("../../src/integration/GaussLegendre.jl")
include("../../src/QuarkDistribution.jl")
include("../../src/QuarkDistribution_Aniso.jl")
include("../../src/relaxtime/OneLoopIntegrals.jl")
include("../../src/relaxtime/OneLoopIntegralsAniso.jl")
include("../../src/relaxtime/PolarizationAniso.jl")

using .GaussLegendre: gauleg
using .PNJLQuarkDistributions_Aniso: correction_cos_theta_coefficient
using .OneLoopIntegralsCorrection: B0_correction
using .PolarizationAniso: polarization_aniso
using .OneLoopIntegrals: A

@testset "ξ<0极化函数问题诊断" begin
    
    # 测试参数
    p = 0.5  # fm⁻¹
    m = 0.25  # fm⁻¹
    μ = 0.12  # fm⁻¹
    T = 0.17  # fm⁻¹
    Φ = 0.15
    Φbar = 0.15
    
    println("\n" * "="^70)
    println("诊断测试1: correction_cos_theta_coefficient 对 ξ 符号的响应")
    println("="^70)
    
    ξ_values = [-0.5, -0.3, -0.1, 0.0, 0.1, 0.3, 0.5]
    
    println("\n参数: p=$p, m=$m, μ=$μ, T=$T")
    println("-"^70)
    println(@sprintf("%-8s | %-20s | %-20s", "ξ", "quark_coeff", "antiquark_coeff"))
    println("-"^70)
    
    for ξ in ξ_values
        quark_coeff = correction_cos_theta_coefficient(:quark, p, m, μ, T, Φ, Φbar, ξ)
        antiquark_coeff = correction_cos_theta_coefficient(:antiquark, p, m, μ, T, Φ, Φbar, ξ)
        
        println(@sprintf("%-8.2f | %-20.12e | %-20.12e", ξ, quark_coeff, antiquark_coeff))
        
        # 检查ξ=0时应该为0
        if abs(ξ) < 1e-10
            @test abs(quark_coeff) < 1e-14
            @test abs(antiquark_coeff) < 1e-14
        end
        
        # 检查线性关系: coeff(ξ) ∝ ξ
        if abs(ξ) > 0.01
            coeff_at_unit_xi = correction_cos_theta_coefficient(:quark, p, m, μ, T, Φ, Φbar, 1.0)
            expected = ξ * coeff_at_unit_xi
            @test isapprox(quark_coeff, expected, rtol=1e-10)
        end
    end
    
    println("\n✓ correction_cos_theta_coefficient 正确响应 ξ 符号（线性关系）")
    
    # ====================================================================
    
    println("\n" * "="^70)
    println("诊断测试2: B0_correction 对 ξ 符号的响应")
    println("="^70)
    
    λ = 0.3  # fm⁻¹
    k = 0.2  # fm⁻¹
    m1 = 0.25  # fm⁻¹
    m2 = 0.30  # fm⁻¹
    μ1 = 0.12  # fm⁻¹
    μ2 = -0.05  # fm⁻¹
    
    println("\n参数: λ=$λ, k=$k, m1=$m1, m2=$m2, μ1=$μ1, μ2=$μ2, T=$T")
    println("-"^70)
    println(@sprintf("%-8s | %-20s | %-20s", "ξ", "B0_real", "B0_imag"))
    println("-"^70)
    
    B0_results = Dict{Float64, Tuple{Float64, Float64}}()
    
    for ξ in ξ_values
        B0_real, B0_imag = B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ)
        B0_results[ξ] = (B0_real, B0_imag)
        
        println(@sprintf("%-8.2f | %-20.12e | %-20.12e", ξ, B0_real, B0_imag))
    end
    
    # 检查ξ=0时是否为基准
    B0_zero = B0_results[0.0]
    println("\n基准值 (ξ=0): Re=$(B0_zero[1]), Im=$(B0_zero[2])")
    
    # 检查正负ξ是否有差异
    println("\n正负ξ对比:")
    println("-"^70)
    for ξ_pos in [0.1, 0.3, 0.5]
        ξ_neg = -ξ_pos
        B0_pos = B0_results[ξ_pos]
        B0_neg = B0_results[ξ_neg]
        
        diff_real = B0_pos[1] - B0_neg[1]
        diff_imag = B0_pos[2] - B0_neg[2]
        
        println(@sprintf("ξ = ±%.1f:", ξ_pos))
        println(@sprintf("  B0(+ξ): Re=%.12e, Im=%.12e", B0_pos[1], B0_pos[2]))
        println(@sprintf("  B0(-ξ): Re=%.12e, Im=%.12e", B0_neg[1], B0_neg[2]))
        println(@sprintf("  差异:    ΔRe=%.12e, ΔIm=%.12e", diff_real, diff_imag))
        
        # 检查是否完全相同（问题症状）
        if abs(diff_real) < 1e-14 && abs(diff_imag) < 1e-14
            println("  ❌ 警告：正负ξ结果完全相同！")
            @test_broken false  # 标记为已知bug
        else
            println("  ✓ 正负ξ有差异")
        end
        println()
    end
    
    # ====================================================================
    
    println("\n" * "="^70)
    println("诊断测试3: polarization_aniso 对 ξ 符号的响应")
    println("="^70)
    
    # 准备A函数值
    nodes_p, weights_p = gauleg(0.0, 20.0, 64)
    A1 = A(m1, μ1, T, Φ, Φbar, nodes_p, weights_p)
    A2 = A(m2, μ2, T, Φ, Φbar, nodes_p, weights_p)
    
    k0 = 0.45  # fm⁻¹
    k_norm = 0.30  # fm⁻¹
    
    println("\n参数: k0=$k0, k=$k_norm, m1=$m1, m2=$m2, μ1=$μ1, μ2=$μ2, T=$T")
    println("     A1=$A1, A2=$A2")
    println("-"^70)
    println(@sprintf("%-8s | %-20s | %-20s | %-20s | %-20s", "ξ", "Π_P_real", "Π_P_imag", "Π_S_real", "Π_S_imag"))
    println("-"^70)
    
    Π_results = Dict{Float64, Dict{Symbol, Tuple{Float64, Float64}}}()
    
    for ξ in ξ_values
        Π_P = polarization_aniso(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
        Π_S = polarization_aniso(:S, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
        
        Π_results[ξ] = Dict(:P => Π_P, :S => Π_S)
        
        println(@sprintf("%-8.2f | %-20.12e | %-20.12e | %-20.12e | %-20.12e", 
                ξ, Π_P[1], Π_P[2], Π_S[1], Π_S[2]))
    end
    
    # 检查正负ξ是否有差异
    println("\n正负ξ对比 (极化函数):")
    println("-"^70)
    for ξ_pos in [0.1, 0.3, 0.5]
        ξ_neg = -ξ_pos
        Π_P_pos = Π_results[ξ_pos][:P]
        Π_P_neg = Π_results[ξ_neg][:P]
        
        diff_real = Π_P_pos[1] - Π_P_neg[1]
        diff_imag = Π_P_pos[2] - Π_P_neg[2]
        
        println(@sprintf("ξ = ±%.1f (P通道):", ξ_pos))
        println(@sprintf("  Π(+ξ): Re=%.12e, Im=%.12e", Π_P_pos[1], Π_P_pos[2]))
        println(@sprintf("  Π(-ξ): Re=%.12e, Im=%.12e", Π_P_neg[1], Π_P_neg[2]))
        println(@sprintf("  差异:   ΔRe=%.12e, ΔIm=%.12e", diff_real, diff_imag))
        
        if abs(diff_real) < 1e-14 && abs(diff_imag) < 1e-14
            println("  ❌ 警告：极化函数对ξ符号不敏感！")
            @test_broken false
        else
            println("  ✓ 极化函数响应ξ符号")
        end
        println()
    end
    
    # ====================================================================
    
    println("\n" * "="^70)
    println("诊断测试4: 检查是否存在abs(ξ)或条件判断")
    println("="^70)
    
    # 手动检查代码中可能的问题点
    println("\n请人工检查以下代码位置:")
    println("1. src/QuarkDistribution_Aniso.jl:74")
    println("   coeff = 0.5 * ξ * (p_inv_fm^2) / E_inv_fm")
    println("   → 检查: 是否误用 abs(ξ)?")
    println()
    println("2. src/relaxtime/OneLoopIntegralsAniso.jl")
    println("   → 检查: B0_correction 各分支中ξ的传递")
    println()
    println("3. src/relaxtime/PolarizationAniso.jl")
    println("   → 检查: polarization_aniso 中是否有 if ξ < 0 的条件")
    println()
end

println("\n" * "="^70)
println("诊断测试完成！")
println("="^70)

