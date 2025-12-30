using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .PNJL.PhaseTransition

# ============================================================================
# 测试数据
# ============================================================================

"""生成 S 形曲线测试数据"""
function s_shape_curve()
    # μ(ρ) 曲线：先增后减再增（S 形）
    rho = Float64[0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    mu = Float64[250.0, 260.0, 270.0, 275.0, 272.0, 268.0, 265.0, 270.0, 280.0, 290.0, 300.0]
    return mu, rho
end

"""生成单调曲线测试数据（无 S 形）"""
function monotonic_curve()
    rho = Float64[0.0, 0.5, 1.0, 1.5, 2.0]
    mu = Float64[250.0, 260.0, 270.0, 280.0, 290.0]
    return mu, rho
end

"""生成合成曲线用于 Maxwell 构造测试"""
function synthetic_maxwell_curve()
    # 更明显的 S 形，便于 Maxwell 构造
    rho = Float64[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
    mu = Float64[250.0, 255.0, 260.0, 265.0, 268.0, 266.0, 262.0, 258.0, 260.0, 265.0, 272.0, 280.0, 288.0, 296.0, 305.0]
    return mu, rho
end

# ============================================================================
# S 形检测测试
# ============================================================================

@testset "detect_s_shape" begin
    @testset "检测 S 形曲线" begin
        mu, rho = s_shape_curve()
        result = detect_s_shape(mu, rho)
        
        @test result.has_s_shape
        @test result.mu_spinodal_hadron !== nothing
        @test result.mu_spinodal_quark !== nothing
        @test result.rho_spinodal_hadron !== nothing
        @test result.rho_spinodal_quark !== nothing
        # hadron 侧是极大值，quark 侧是极小值
        @test result.mu_spinodal_hadron > result.mu_spinodal_quark
        @test result.rho_spinodal_hadron < result.rho_spinodal_quark
    end
    
    @testset "单调曲线无 S 形" begin
        mu, rho = monotonic_curve()
        result = detect_s_shape(mu, rho)
        
        @test !result.has_s_shape
        @test result.mu_spinodal_hadron === nothing
        @test result.mu_spinodal_quark === nothing
    end
    
    @testset "数据点不足" begin
        mu = Float64[250.0, 260.0]
        rho = Float64[0.0, 1.0]
        result = detect_s_shape(mu, rho; min_points=5)
        
        @test !result.has_s_shape
    end
end

# ============================================================================
# Maxwell 构造测试
# ============================================================================

@testset "maxwell_construction" begin
    @testset "基本 Maxwell 构造" begin
        mu, rho = synthetic_maxwell_curve()
        result = maxwell_construction(mu, rho; min_samples=8)
        
        @test result.converged
        @test result.mu_transition !== nothing
        @test result.rho_hadron !== nothing
        @test result.rho_quark !== nothing
        @test result.rho_hadron < result.rho_quark
        @test result.area_residual !== nothing
        @test result.area_residual < 0.1  # 面积残差应该很小
    end
    
    @testset "使用 spinodal_hint" begin
        mu, rho = synthetic_maxwell_curve()
        hint = detect_s_shape(mu, rho)
        result = maxwell_construction(mu, rho; spinodal_hint=hint, min_samples=8)
        
        @test result.converged
    end
    
    @testset "单调曲线无法构造" begin
        mu, rho = monotonic_curve()
        result = maxwell_construction(mu, rho; min_samples=3)
        
        @test !result.converged
        @test get(result.details, :reason, "") == "no_s_shape"
    end
    
    @testset "数据点不足" begin
        mu = Float64[250.0, 260.0, 270.0]
        rho = Float64[0.0, 0.5, 1.0]
        result = maxwell_construction(mu, rho; min_samples=12)
        
        @test !result.converged
        @test get(result.details, :reason, "") == "insufficient_points"
    end
end

# ============================================================================
# 工具函数测试
# ============================================================================

@testset "group_curves_by_temperature" begin
    rows = [
        Dict("T_MeV" => "70", "xi" => "0.0", "rho" => "0.1", "mu_avg_MeV" => "250.0"),
        Dict("T_MeV" => "70", "xi" => "0.0", "rho" => "0.2", "mu_avg_MeV" => "260.0"),
        Dict("T_MeV" => "80", "xi" => "0.0", "rho" => "0.1", "mu_avg_MeV" => "255.0"),
        Dict("T_MeV" => "70", "xi" => "0.5", "rho" => "0.1", "mu_avg_MeV" => "245.0"),  # 不同 xi
    ]
    
    grouped = group_curves_by_temperature(rows; xi=0.0)
    
    @test haskey(grouped, 70.0)
    @test haskey(grouped, 80.0)
    @test length(grouped[70.0]) == 2  # xi=0.5 的被过滤
    @test length(grouped[80.0]) == 1
end
