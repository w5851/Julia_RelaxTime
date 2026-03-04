# LegacyPNJLModel.jl 单元测试
#
# 测试内容：
# 1. LegacyPNJLModel 构造
# 2. 类型继承
# 3. calculate_mass_vec 委托

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# ============================================================================

@testset "LegacyPNJLModel" begin

    @testset "构造" begin
        m = Models.LegacyPNJLModel()
        @test m isa Models.AbstractPNJLModel
        @test m isa Models.LegacyPNJLModel
    end

    @testset "create_model(:PNJL) 不返回 Legacy" begin
        m = Models.create_model(:PNJL)
        @test !(m isa Models.LegacyPNJLModel)
        @test m isa Models.PNJLModel
    end

    @testset "calculate_mass_vec 委托" begin
        m = Models.LegacyPNJLModel()
        φ = SVector{3}(-1.84, -1.84, -2.23)
        masses = Base.invokelatest(Models.calculate_mass_vec, m, φ)
        @test length(masses) == 3
        @test all(isfinite.(masses))
        @test all(masses .> 0)
    end
end
