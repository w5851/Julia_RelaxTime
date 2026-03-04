# Models abstract_model.jl 单元测试
#
# 测试内容：
# 1. 抽象类型层级
# 2. 接口默认抛出 MethodError
# 3. gap_state_dim / gap_residual 契约

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

# ============================================================================
# 抽象类型层级
# ============================================================================

@testset "Abstract model types" begin

    @testset "类型层级" begin
        @test Models.AbstractNJLModel <: Models.AbstractQCDModel
        @test Models.AbstractPNJLModel <: Models.AbstractNJLModel
        @test Models.AbstractPNJLModel <: Models.AbstractQCDModel
    end

    @testset "NJLModel 是 AbstractNJLModel" begin
        m = Models.create_model(:NJL)
        @test m isa Models.AbstractNJLModel
        @test m isa Models.AbstractQCDModel
    end

    @testset "PNJLModel 是 AbstractPNJLModel" begin
        Models.pnjl_module()
        m = Models.create_model(:PNJL)
        @test m isa Models.AbstractPNJLModel
        @test m isa Models.AbstractNJLModel
        @test m isa Models.AbstractQCDModel
    end
end

# ============================================================================
# 默认抛出 MethodError
# ============================================================================

@testset "Default interface throws MethodError" begin
    # 创建一个 bare 子类型来测试默认行为
    struct _TestDummyModel <: Models.AbstractQCDModel end

    dummy = _TestDummyModel()
    T = 0.5
    μ = SVector{3}(0.0, 0.0, 0.0)
    φ = SVector{3}(0.1, 0.2, 0.3)

    @test_throws MethodError Models.solve_gap(dummy, T, μ)
    @test_throws MethodError Models.number_densities(dummy, φ, T, μ)
    @test_throws MethodError Models.calculate_mass_vec(dummy, φ)
    @test_throws MethodError Models.calculate_chiral(dummy, φ)
end

# ============================================================================
# 接口函数可用性
# ============================================================================

@testset "Interface functions defined" begin
    @test isdefined(Models, :solve_gap)
    @test isdefined(Models, :number_densities)
    @test isdefined(Models, :calculate_mass_vec)
    @test isdefined(Models, :calculate_chiral)
    @test isdefined(Models, :gap_state_dim)
    @test isdefined(Models, :gap_residual)
end
