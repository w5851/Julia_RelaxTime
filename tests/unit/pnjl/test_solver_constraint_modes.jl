# PNJL Solver ConstraintModes 单元测试
#
# 测试内容：
# 1. 约束模式定义
# 2. 状态维度
# 3. 参数维度
# 4. 描述字符串

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

# 避免跨测试文件导出冲突：通过模块前缀访问。
P = PNJL

# ============================================================================
# 模式构造测试
# ============================================================================

@testset "ConstraintModes construction" begin
    @testset "FixedMu" begin
        mode = P.FixedMu()
        @test mode isa P.ConstraintMode
        @test mode isa P.FixedMu
    end
    
    @testset "FixedRho" begin
        mode = P.FixedRho(1.0)
        @test mode isa P.ConstraintMode
        @test mode isa P.FixedRho
        @test mode.rho_target == 1.0
    end
    
    @testset "FixedEntropy" begin
        mode = P.FixedEntropy(0.5)
        @test mode isa P.ConstraintMode
        @test mode isa P.FixedEntropy
        @test mode.s_target == 0.5
    end
    
    @testset "FixedSigma" begin
        mode = P.FixedSigma(10.0)
        @test mode isa P.ConstraintMode
        @test mode isa P.FixedSigma
        @test mode.sigma_target == 10.0
    end
end

# ============================================================================
# 维度查询测试
# ============================================================================

@testset "state_dim" begin
    @test P.state_dim(P.FixedMu()) == 5
    @test P.state_dim(P.FixedRho(1.0)) == 8
    @test P.state_dim(P.FixedEntropy(0.5)) == 8
    @test P.state_dim(P.FixedSigma(10.0)) == 8
end

@testset "param_dim" begin
    @test P.param_dim(P.FixedMu()) == 2  # [T, μ]
    @test P.param_dim(P.FixedRho(1.0)) == 1  # [T]
    @test P.param_dim(P.FixedEntropy(0.5)) == 1  # [T]
    @test P.param_dim(P.FixedSigma(10.0)) == 1  # [T]
end

# ============================================================================
# 描述字符串测试
# ============================================================================

@testset "constraint_description" begin
    @test P.constraint_description(P.FixedMu()) == "Fixed chemical potential μ"
    @test P.constraint_description(P.FixedRho(1.5)) == "Fixed baryon density ρ/ρ₀ = 1.5"
    @test P.constraint_description(P.FixedEntropy(0.3)) == "Fixed entropy density s = 0.3 fm⁻³"
    @test P.constraint_description(P.FixedSigma(8.0)) == "Fixed specific entropy σ = s/n_B = 8.0"
end

# ============================================================================
# 显示方法测试
# ============================================================================

@testset "show methods" begin
    @test sprint(show, P.FixedMu()) == "FixedMu()"
    @test sprint(show, P.FixedRho(1.0)) == "FixedRho(ρ/ρ₀=1.0)"
    @test sprint(show, P.FixedEntropy(0.5)) == "FixedEntropy(s=0.5)"
    @test sprint(show, P.FixedSigma(10.0)) == "FixedSigma(σ=10.0)"
end
