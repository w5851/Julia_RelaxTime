# 测试 MomentumMapping 模块

using Test
using LinearAlgebra

# 添加模块路径
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src", "simulation"))

# 导入模块
include("../src/simulation/FrameTransformations.jl")
include("../src/simulation/EllipsoidCalculation.jl")
include("../src/simulation/MomentumMapping.jl")

using .MomentumMapping

@testset "MomentumMapping Tests" begin
    
    @testset "Basic scattering calculation" begin
        # 测试默认参数的散射计算
        p1 = [0.5, 0.0, 1.8]
        p2 = [-0.5, 0.0, -1.8]
        m1 = m2 = m3 = m4 = 1.52
        
        result = calculate_outgoing_momenta(p1, p2, m1, m2, m3, m4)
        
        # 基本检查
        @test !isnothing(result)
        @test result.s > 0
        @test result.p_star > 0
        @test result.gamma >= 1.0
        @test 0 <= norm(result.beta) < 1.0
        
        println("✓ Basic calculation successful:")
        println("  √s = ", sqrt(result.s))
        println("  p* = ", result.p_star)
    end
    
    @testset "Physics validation" begin
        p1 = [0.5, 0.0, 1.8]
        p2 = [-0.5, 0.0, -1.8]
        m = 1.52
        
        result = calculate_outgoing_momenta(p1, p2, m, m, m, m, π/4, π/6)
        is_valid, checks = validate_kinematics(result, m, m, m, m)
        
        println("\n✓ Physics validation:")
        println("  Energy conservation: ", checks["energy_conservation"])
        println("  Momentum conservation: ", checks["momentum_conservation"])
        println("  On-shell p3: ", checks["on_shell_p3"])
        println("  On-shell p4: ", checks["on_shell_p4"])
        
        # 应该满足守恒律
        @test is_valid
        @test checks["energy_conservation"][1] < 1e-9
        @test checks["momentum_conservation"][1] < 1e-9
    end
    
    @testset "Different angles" begin
        p1 = [0.5, 0.0, 1.8]
        p2 = [-0.5, 0.0, -1.8]
        m = 1.52
        
        # 测试多个角度
        angles = [
            (0.0, 0.0),        # 前向散射
            (π, 0.0),          # 后向散射
            (π/2, 0.0),        # 90度散射
            (π/4, π/4),        # 对角散射
        ]
        
        for (theta, phi) in angles
            result = calculate_outgoing_momenta(p1, p2, m, m, m, m, theta, phi)
            is_valid, _ = validate_kinematics(result, m, m, m, m)
            
            @test is_valid
            println("  θ=$(theta/π)π, φ=$(phi/π)π: ✓")
        end
    end
    
    @testset "Ellipsoid properties" begin
        p1 = [0.5, 0.0, 1.8]
        p2 = [-0.5, 0.0, -1.8]
        m = 1.52
        
        result = calculate_outgoing_momenta(p1, p2, m, m, m, m)
        ellipsoid = result.ellipsoid
        
        # 椭球参数应该合理
        @test length(ellipsoid.center) == 3
        @test size(ellipsoid.axes_directions) == (3, 3)
        @test length(ellipsoid.half_lengths) == 3
        @test all(ellipsoid.half_lengths .> 0)
        
        # 主轴方向应该正交
        axes = ellipsoid.axes_directions
        @test isapprox(dot(axes[:, 1], axes[:, 2]), 0.0, atol=1e-9)
        @test isapprox(dot(axes[:, 2], axes[:, 3]), 0.0, atol=1e-9)
        @test isapprox(dot(axes[:, 3], axes[:, 1]), 0.0, atol=1e-9)
        
        println("✓ Ellipsoid properties verified:")
        println("  Center: ", ellipsoid.center)
        println("  Half-lengths: ", ellipsoid.half_lengths)
    end
    
    @testset "Mandelstam variables" begin
        p1 = [0.5, 0.0, 1.8]
        p2 = [-0.5, 0.0, -1.8]
        m = 1.52
        
        result = calculate_outgoing_momenta(p1, p2, m, m, m, m, π/4, π/6)
        
        t = calculate_mandelstam_t(result)
        u = calculate_mandelstam_u(result)
        
        # Mandelstam关系: s + t + u = Σmᵢ²
        sum_masses_sq = m^2 + m^2 + m^2 + m^2
        @test isapprox(result.s + t + u, sum_masses_sq, atol=1e-8)
        
        println("✓ Mandelstam variables:")
        println("  s = ", result.s)
        println("  t = ", t)
        println("  u = ", u)
        println("  s+t+u = ", result.s + t + u, " (expected: ", sum_masses_sq, ")")
    end
    
    @testset "Symmetric scattering" begin
        # 完全对称情况：质心系静止
        p1 = [0.0, 0.0, 2.0]
        p2 = [0.0, 0.0, -2.0]
        m = 1.0
        
        result = calculate_outgoing_momenta(p1, p2, m, m, m, m, π/2, 0.0)
        
        # β应该接近0
        @test norm(result.beta) < 1e-10
        @test result.gamma ≈ 1.0 atol=1e-10
        
        # 椭球中心应该在原点
        @test norm(result.ellipsoid.center) < 1e-10
        
        println("✓ Symmetric scattering verified")
    end
    
    @testset "Threshold check" begin
        # 测试阈值附近的情况
        p1 = [0.0, 0.0, 0.1]  # 很小的动量
        p2 = [0.0, 0.0, -0.1]
        m = 1.0
        
        # 这应该会抛出阈值错误
        @test_throws ErrorException calculate_outgoing_momenta(p1, p2, m, m, m, m)
        
        println("✓ Threshold check working")
    end
end

println("\n" * "="^60)
println("All MomentumMapping tests passed! ✅")
println("="^60)
