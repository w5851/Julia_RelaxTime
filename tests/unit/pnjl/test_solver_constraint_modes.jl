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

using .PNJL.ConstraintModes

# ============================================================================
# 模式构造测试
# ============================================================================

@testset "ConstraintModes construction" begin
    @testset "FixedMu" begin
        mode = FixedMu()
        @test mode isa ConstraintMode
        @test mode isa FixedMu
    end
    
    @testset "FixedRho" begin
        mode = FixedRho(1.0)
        @test mode isa ConstraintMode
        @test mode isa FixedRho
        @test mode.rho_target == 1.0
    end
    
    @testset "FixedEntropy" begin
        mode = FixedEntropy(0.5)
        @test mode isa ConstraintMode
        @test mode isa FixedEntropy
        @test mode.s_target == 0.5
    end
    
    @testset "FixedSigma" begin
        mode = FixedSigma(10.0)
        @test mode isa ConstraintMode
        @test mode isa FixedSigma
        @test mode.sigma_target == 10.0
    end
end

# ============================================================================
# 维度查询测试
# ============================================================================

@testset "state_dim" begin
    @test state_dim(FixedMu()) == 5
    @test state_dim(FixedRho(1.0)) == 8
    @test state_dim(FixedEntropy(0.5)) == 8
    @test state_dim(FixedSigma(10.0)) == 8
end

@testset "param_dim" begin
    @test param_dim(FixedMu()) == 2  # [T, μ]
    @test param_dim(FixedRho(1.0)) == 1  # [T]
    @test param_dim(FixedEntropy(0.5)) == 1  # [T]
    @test param_dim(FixedSigma(10.0)) == 1  # [T]
end

# ============================================================================
# 描述字符串测试
# ============================================================================

@testset "constraint_description" begin
    @test constraint_description(FixedMu()) == "Fixed chemical potential μ"
    @test constraint_description(FixedRho(1.5)) == "Fixed baryon density ρ/ρ₀ = 1.5"
    @test constraint_description(FixedEntropy(0.3)) == "Fixed entropy density s = 0.3 fm⁻³"
    @test constraint_description(FixedSigma(8.0)) == "Fixed specific entropy σ = s/n_B = 8.0"
end

# ============================================================================
# 显示方法测试
# ============================================================================

@testset "show methods" begin
    @test sprint(show, FixedMu()) == "FixedMu()"
    @test sprint(show, FixedRho(1.0)) == "FixedRho(ρ/ρ₀=1.0)"
    @test sprint(show, FixedEntropy(0.5)) == "FixedEntropy(s=0.5)"
    @test sprint(show, FixedSigma(10.0)) == "FixedSigma(σ=10.0)"
end
