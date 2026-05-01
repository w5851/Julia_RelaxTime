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

    @testset "meson_density_workflow_module" begin
        mod = Models.meson_density_workflow_module()
        @test mod === Models.MesonDensityWorkflow
    end

    @testset "magnetic_thermodynamics_module" begin
        mod = Models.magnetic_thermodynamics_module()
        @test mod === Models.MagneticThermodynamics
    end

    @testset "积分网格契约 accessor" begin
        @test isdefined(Models, :default_momentum_count)
        @test isdefined(Models, :default_theta_count)
        @test isdefined(Models, :default_momentum_nodes)
        @test isdefined(Models, :default_momentum_weights)

        @test Models.default_momentum_count() == Models.PNJLCore.DEFAULT_MOMENTUM_COUNT
        @test Models.default_theta_count() == Models.PNJLCore.DEFAULT_THETA_COUNT
        @test Models.default_momentum_nodes() === Models.PNJLIntegrals.THERMAL_DEFAULT_NODES
        @test Models.default_momentum_weights() === Models.PNJLIntegrals.THERMAL_DEFAULT_WEIGHTS
    end

    @testset "薄转发接口存在" begin
        @test isdefined(Models, :run_tmu_scan)
        @test isdefined(Models, :run_trho_scan)
        @test isdefined(Models, :build_default_rho_grid)
        @test isdefined(Models, :default_scan_numeric_options)
        @test isdefined(Models, :solve_pnjl_point)
        @test isdefined(Models, :solve_gap_and_transport)
        @test isdefined(Models, :solve_transport_from_equilibrium)
        @test isdefined(Models, :solve_gap_and_meson_point)
        @test isdefined(Models, :solve_meson_density_from_meson_point)
        @test isdefined(Models, :solve_gap_and_meson_density_point)
        @test isdefined(Models, :run_phase_pipeline)
        @test isdefined(Models, :run_production_phase_pipeline)
        @test isdefined(Models, :find_cep)
        @test isdefined(Models, :build_phase_artifacts)
        @test isdefined(Models, :resolve_phase_output_target)
        @test isdefined(Models, :promote_phase_artifacts)
        @test isdefined(Models, :auto_phase_hint)
    end

    @testset "同步/异步入口默认数值参数对齐" begin
        opts = Models.default_scan_numeric_options()
        @test opts.p_num == 24
        @test opts.t_num == 8
    end

    @testset "solve_pnjl_point 参数校验" begin
        @test_throws ArgumentError Models.solve_pnjl_point(T_mev=150.0)
        @test_throws ArgumentError Models.solve_pnjl_point(T_mev=NaN, mu_mev=0.0)
    end
end
