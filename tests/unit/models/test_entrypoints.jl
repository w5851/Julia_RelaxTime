# entrypoints.jl 单元测试
#
# 测试内容：
# 1. 薄转发函数接口存在
# 2. 模块访问器返回正确模块

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# ============================================================================

@testset "entrypoints" begin

    @testset "pnjl_module 返回 Models" begin
        mod = Models.pnjl_module()
        @test mod === Models
    end

    @testset "transport_workflow_module" begin
        mod = Models.transport_workflow_module()
        @test mod === Models.TransportWorkflow
    end

    @testset "meson_workflow_module" begin
        mod = Models.meson_workflow_module()
        @test mod === Models.MesonMassWorkflow
    end

    @testset "薄转发接口存在" begin
        @test isdefined(Models, :run_tmu_scan)
        @test isdefined(Models, :run_trho_scan)
        @test isdefined(Models, :build_default_rho_grid)
        @test isdefined(Models, :solve_gap_and_transport)
        @test isdefined(Models, :solve_transport_from_equilibrium)
        @test isdefined(Models, :solve_gap_and_meson_point)
        @test isdefined(Models, :run_phase_pipeline)
        @test isdefined(Models, :find_cep)
        @test isdefined(Models, :build_phase_artifacts)
        @test isdefined(Models, :resolve_phase_output_target)
        @test isdefined(Models, :promote_phase_artifacts)
    end
end
