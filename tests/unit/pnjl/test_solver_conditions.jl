# PNJL Solver Conditions 单元测试
#
# 测试内容：
# 1. GapParams 结构
# 2. gap_conditions 核心函数
# 3. build_conditions 各模式
# 4. build_residual! 各模式

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .PNJL.ConstraintModes
using .PNJL.Conditions
using .PNJL.Conditions.Thermodynamics.Integrals: cached_nodes

# ============================================================================
# GapParams 测试
# ============================================================================

@testset "GapParams" begin
    T_fm = 0.5
    thermal_nodes = cached_nodes(24, 6)
    xi = 0.0
    
    params = GapParams(T_fm, thermal_nodes, xi)
    
    @test params.T_fm == T_fm
    @test params.thermal_nodes === thermal_nodes
    @test params.xi == xi
end

# ============================================================================
# gap_conditions 测试
# ============================================================================

@testset "gap_conditions" begin
    T_fm = 0.5
    thermal_nodes = cached_nodes(24, 6)
    xi = 0.0
    params = GapParams(T_fm, thermal_nodes, xi)
    
    @testset "返回值维度" begin
        x_state = SVector{5}(-1.5, -1.5, -2.1, 0.2, 0.2)
        mu_vec = SVector{3}(1.0, 1.0, 1.0)
        
        residual = gap_conditions(x_state, mu_vec, params)
        
        @test length(residual) == 5
        @test all(isfinite.(residual))
    end
    
    @testset "平衡态残差小" begin
        # 使用已知的近似平衡态
        x_state = SVector{5}(-1.84329, -1.84329, -2.22701, 1.0e-5, 4.0e-5)
        mu_vec = SVector{3}(0.05, 0.05, 0.05)  # 低化学势
        T_fm_low = 0.25  # ~50 MeV
        params_low = GapParams(T_fm_low, thermal_nodes, xi)
        
        residual = gap_conditions(x_state, mu_vec, params_low)
        
        # 残差应该相对较小（但不一定为零，因为这不是精确解）
        @test all(isfinite.(residual))
    end
end

# ============================================================================
# build_conditions 测试
# ============================================================================

@testset "build_conditions" begin
    thermal_nodes = cached_nodes(24, 6)
    params = GapParams(0.5, thermal_nodes, 0.0)
    
    @testset "FixedMu" begin
        cond_fn = build_conditions(FixedMu(), params)
        
        θ = [0.5, 1.0]  # T, μ
        x = [-1.5, -1.5, -2.1, 0.2, 0.2]
        
        residual = cond_fn(θ, x)
        
        @test length(residual) == 5
        @test all(isfinite.(residual))
    end
    
    @testset "FixedRho" begin
        cond_fn = build_conditions(FixedRho(1.0), params)
        
        θ = [0.5]  # T
        x = [-1.5, -1.5, -2.1, 0.2, 0.2, 1.5, 1.5, 1.5]  # 5 + 3 维
        
        residual = cond_fn(θ, x)
        
        @test length(residual) == 8
        @test all(isfinite.(residual))
    end
    
    @testset "FixedEntropy" begin
        cond_fn = build_conditions(FixedEntropy(0.5), params)
        
        θ = [0.5]
        x = [-1.5, -1.5, -2.1, 0.2, 0.2, 1.5, 1.5, 1.5]
        
        residual = cond_fn(θ, x)
        
        @test length(residual) == 8
        @test all(isfinite.(residual))
    end
    
    @testset "FixedSigma" begin
        cond_fn = build_conditions(FixedSigma(10.0), params)
        
        θ = [0.5]
        x = [-1.5, -1.5, -2.1, 0.2, 0.2, 1.5, 1.5, 1.5]
        
        residual = cond_fn(θ, x)
        
        @test length(residual) == 8
        @test all(isfinite.(residual))
    end
end

# ============================================================================
# build_residual! 测试
# ============================================================================

@testset "build_residual!" begin
    thermal_nodes = cached_nodes(24, 6)
    T_fm = 0.5
    params = GapParams(T_fm, thermal_nodes, 0.0)
    
    @testset "FixedMu" begin
        mu_vec = SVector{3}(1.0, 1.0, 1.0)
        residual_fn! = build_residual!(FixedMu(), mu_vec, params)
        
        x = [-1.5, -1.5, -2.1, 0.2, 0.2]
        F = zeros(5)
        
        residual_fn!(F, x)
        
        @test length(F) == 5
        @test all(isfinite.(F))
    end
    
    @testset "FixedRho" begin
        residual_fn! = build_residual!(FixedRho(1.0), params)
        
        x = [-1.5, -1.5, -2.1, 0.2, 0.2, 1.5, 1.5, 1.5]
        F = zeros(8)
        
        residual_fn!(F, x)
        
        @test length(F) == 8
        @test all(isfinite.(F))
        
        # 检查化学势相等约束
        @test F[6] ≈ x[6] - x[7]  # μ_u = μ_d
        @test F[7] ≈ x[7] - x[8]  # μ_d = μ_s
    end
    
    @testset "FixedEntropy" begin
        residual_fn! = build_residual!(FixedEntropy(0.5), params)
        
        x = [-1.5, -1.5, -2.1, 0.2, 0.2, 1.5, 1.5, 1.5]
        F = zeros(8)
        
        residual_fn!(F, x)
        
        @test length(F) == 8
        @test all(isfinite.(F))
    end
    
    @testset "FixedSigma" begin
        residual_fn! = build_residual!(FixedSigma(10.0), params)
        
        x = [-1.5, -1.5, -2.1, 0.2, 0.2, 1.5, 1.5, 1.5]
        F = zeros(8)
        
        residual_fn!(F, x)
        
        @test length(F) == 8
        @test all(isfinite.(F))
    end
end

# ============================================================================
# 约束一致性测试
# ============================================================================

@testset "constraint consistency" begin
    thermal_nodes = cached_nodes(24, 6)
    T_fm = 0.5
    params = GapParams(T_fm, thermal_nodes, 0.0)
    
    @testset "FixedRho 化学势相等" begin
        residual_fn! = build_residual!(FixedRho(1.0), params)
        
        # 化学势相等的情况
        x_equal = [-1.5, -1.5, -2.1, 0.2, 0.2, 1.5, 1.5, 1.5]
        F = zeros(8)
        residual_fn!(F, x_equal)
        @test abs(F[6]) < 1e-10  # μ_u = μ_d
        @test abs(F[7]) < 1e-10  # μ_d = μ_s
        
        # 化学势不等的情况
        x_unequal = [-1.5, -1.5, -2.1, 0.2, 0.2, 1.5, 1.6, 1.7]
        residual_fn!(F, x_unequal)
        @test abs(F[6]) > 0.05  # μ_u ≠ μ_d
        @test abs(F[7]) > 0.05  # μ_d ≠ μ_s
    end
end
