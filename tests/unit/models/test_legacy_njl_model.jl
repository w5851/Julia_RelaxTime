# LegacyNJLModel.jl 单元测试
#
# 测试内容：
# 1. LegacyNJLModel 构造
# 2. 类型继承
# 3. 委托到内部 NJLModel

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

# ============================================================================

@testset "LegacyNJLModel" begin

    @testset "构造" begin
        m = Models.LegacyNJLModel()
        @test m isa Models.AbstractNJLModel
        @test m isa Models.LegacyNJLModel
        @test m.inner isa Models.NJLModel
    end

    @testset "calculate_mass_vec 委托" begin
        m = Models.LegacyNJLModel()
        φ = SVector{3}(-1.84, -1.84, -2.23)
        masses = Models.calculate_mass_vec(m, φ)
        # 应与 NJLModel 一致
        masses_njl = Models.calculate_mass_vec(m.inner, φ)
        @test all(masses .≈ masses_njl)
    end

    @testset "vacuum_contribution 委托" begin
        m = Models.LegacyNJLModel()
        φ = SVector{3}(-1.84, -1.84, -2.23)
        masses = Models.calculate_mass_vec(m, φ)
        vac = Models.vacuum_contribution(m, masses)
        vac_njl = Models.vacuum_contribution(m.inner, masses)
        @test vac ≈ vac_njl
    end
end
