# MesonMassWorkflow.jl 单元测试
#
# 测试内容：
# 1. 模块加载
# 2. DEFAULT_MESONS 常量
# 3. solve_gap_and_meson_point / build_equilibrium_params 接口存在

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

const _MW = Models.MesonMassWorkflow

# ============================================================================

@testset "MesonMassWorkflow" begin

    @testset "DEFAULT_MESONS 常量" begin
        @test isdefined(_MW, :DEFAULT_MESONS)
        mesons = _MW.DEFAULT_MESONS
        @test mesons isa Tuple
        @test :pi in mesons
        @test :K in mesons
        @test :sigma_pi in mesons
        @test length(mesons) >= 6
    end

    @testset "solve_gap_and_meson_point 接口存在" begin
        @test isdefined(_MW, :solve_gap_and_meson_point)
        @test _MW.solve_gap_and_meson_point isa Function
    end

    @testset "build_equilibrium_params 接口存在" begin
        @test isdefined(_MW, :build_equilibrium_params)
        @test _MW.build_equilibrium_params isa Function
    end
end
