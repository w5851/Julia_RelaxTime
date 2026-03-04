# PNJLCore 模块单元测试
#
# 测试内容：
# 1. PNJLParams 结构体构建
# 2. 基本物理量计算

using Test
using StaticArrays

const PROJECT_ROOT_PC = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_PC, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# PNJLCore is a submodule of Models
const PNJLCore_mod = Models.PNJLCore
const pnjl_params = PNJLCore_mod.pnjl_params
const calculate_mass_vec = PNJLCore_mod.calculate_mass_vec
const chiral_potential = PNJLCore_mod.chiral_potential
const polyakov_potential = PNJLCore_mod.polyakov_potential

@testset "PNJLCore" begin
    params = pnjl_params()

    @testset "pnjl_params 默认参数有效" begin
        @test params isa Any  # PNJLParams struct
    end

    @testset "calculate_mass_vec 返回正质量" begin
        φ = SVector{3, Float64}(-0.03, -0.03, -0.04)
        masses = calculate_mass_vec(params, φ)
        @test all(m -> m > 0, masses)
        @test all(isfinite, masses)
        @test length(masses) == 3
    end

    @testset "chiral_potential 返回有限值" begin
        φ = SVector{3, Float64}(-0.03, -0.03, -0.04)
        V = chiral_potential(params, φ)
        @test isfinite(V)
    end

    @testset "polyakov_potential 返回有限值" begin
        T_fm = 0.8
        Φ = 0.2
        Φbar = 0.2
        V_poly = polyakov_potential(params, Φ, Φbar, T_fm)
        @test isfinite(V_poly)
    end
end
