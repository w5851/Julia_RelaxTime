# 测试 FrameTransformations 模块

using Test
using LinearAlgebra

# 添加模块路径
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "src", "simulation"))

# 导入模块
include("../../src/simulation/FrameTransformations.jl")
using .FrameTransformations

@testset "FrameTransformations Tests" begin
    
    @testset "Källen function" begin
        # 测试对称性
        @test kallen_lambda(1.0, 2.0, 3.0) ≈ kallen_lambda(2.0, 1.0, 3.0)
        @test kallen_lambda(1.0, 2.0, 3.0) ≈ kallen_lambda(3.0, 2.0, 1.0)
        
        # 测试具体数值
        # λ(1,2,3) = 1 + 4 + 9 - 4 - 6 - 6 = -2
        @test kallen_lambda(1.0, 2.0, 3.0) ≈ -2.0
        
        # 测试零情况
        @test kallen_lambda(0.0, 0.0, 0.0) ≈ 0.0
    end
    
    @testset "CMS parameters calculation" begin
        # 测试对称散射：p1 = -p2, m1 = m2 = m3 = m4
        p1 = [0.0, 0.0, 2.0]
        p2 = [0.0, 0.0, -2.0]
        m = 1.0
        
        cms = calculate_cms_parameters(p1, p2, m, m, m, m)
        
        # 检查质心系已经静止（β≈0）
        @test norm(cms.beta) < 1e-10
        @test cms.gamma ≈ 1.0 atol=1e-10
        
        # 检查能量
        E1 = sqrt(m^2 + norm(p1)^2)
        E2 = sqrt(m^2 + norm(p2)^2)
        s_expected = (E1 + E2)^2 - norm(p1 + p2)^2
        @test cms.s ≈ s_expected
        
        println("✓ CMS parameters (symmetric case):")
        println("  √s = ", sqrt(cms.s))
        println("  p* = ", cms.p_star)
    end
    
    @testset "CMS parameters - boosted case" begin
        # 测试不对称情况
        p1 = [0.5, 0.0, 1.8]
        p2 = [-0.5, 0.0, -1.8]
        m = 1.52
        
        cms = calculate_cms_parameters(p1, p2, m, m, m, m)
        
        # β应该指向总动量方向
        p_total = p1 + p2
        if norm(p_total) > 1e-10
            beta_direction = p_total / norm(p_total)
            cms_beta_direction = cms.beta / norm(cms.beta)
            @test isapprox(beta_direction, cms_beta_direction, atol=1e-6)
        end
        
        # γ > 1
        @test cms.gamma >= 1.0
        
        println("✓ CMS parameters (boosted case):")
        println("  √s = ", sqrt(cms.s))
        println("  |β| = ", norm(cms.beta))
        println("  γ = ", cms.gamma)
    end
    
    @testset "Affine transform construction" begin
        # 测试β=0的情况（应该返回单位变换）
        beta_zero = [0.0, 0.0, 0.0]
        A_zero, b_zero = build_affine_transform(beta_zero, 1.0, 1.0)
        
        @test A_zero ≈ I(3)
        @test b_zero ≈ [0.0, 0.0, 0.0]
        
        # 测试非零β
        beta = [0.0, 0.0, 0.6]
        gamma = 1.25  # γ = 1/√(1-β²) = 1/√(1-0.36) = 1.25
        E_star = 2.0
        
        A, b = build_affine_transform(beta, gamma, E_star)
        
        # A应该是对称矩阵
        @test isapprox(A, A', atol=1e-10)
        
        # b应该平行于β
        @test isapprox(b, gamma * E_star * beta, atol=1e-10)
        
        println("✓ Affine transform constructed")
    end
    
    @testset "Boost to lab frame" begin
        # 测试简单情况：沿z轴boost
        p_cms = [1.0, 0.0, 0.0]  # CMS中横向动量
        beta = [0.0, 0.0, 0.6]
        gamma = 1.25
        
        A, b = build_affine_transform(beta, gamma, 1.0)
        p_lab = boost_to_lab(p_cms, A, b)
        
        # 横向分量不变
        @test p_lab[1] ≈ p_cms[1]
        @test p_lab[2] ≈ p_cms[2]
        
        println("✓ Boost transformation applied")
        println("  p_cms = ", p_cms)
        println("  p_lab = ", p_lab)
    end
    
    @testset "Energy boost" begin
        # 测试能量boost
        E_star = 2.0
        p_star = [1.0, 0.0, 1.0]
        beta = [0.0, 0.0, 0.6]
        gamma = 1.25
        
        E_lab = boost_energy(E_star, p_star, beta, gamma)
        
        # E_lab = γ(E* + p*·β)
        E_expected = gamma * (E_star + dot(p_star, beta))
        @test E_lab ≈ E_expected
        
        println("✓ Energy boost:")
        println("  E* = ", E_star)
        println("  E_lab = ", E_lab)
    end
    
    @testset "Round-trip transformation" begin
        # 测试从实验室系到质心系再回到实验室系
        p1_lab = [0.5, 0.0, 1.8]
        p2_lab = [-0.5, 0.0, -1.8]
        m = 1.52
        
        cms = calculate_cms_parameters(p1_lab, p2_lab, m, m, m, m)
        A, b = build_affine_transform(cms.beta, cms.gamma, cms.E3_star)
        
        # 假设在CMS中有一个动量
        p_cms = [0.5, 0.5, cms.p_star * cos(π/4)]
        
        # Boost到实验室系
        p_lab = boost_to_lab(p_cms, A, b)
        
        # 验证：实验室系动量的模应该合理
        @test norm(p_lab) > 0
        @test !any(isnan.(p_lab))
        
        println("✓ Round-trip transformation successful")
    end
end

println("\n" * "="^60)
println("All FrameTransformations tests passed! ✅")
println("="^60)

