# PhasePipeline.jl 单元测试
#
# 测试内容：
# 1. run_phase_pipeline 接口存在
# 2. PhasePipelineResult 结构字段

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

# ============================================================================

@testset "PhasePipeline" begin

    @testset "run_phase_pipeline 接口存在" begin
        @test isdefined(Models, :run_phase_pipeline)
        @test Models.run_phase_pipeline isa Function
    end

    @testset "PhasePipelineResult 字段完整" begin
        @test isdefined(Models, :PhasePipelineResult)
        # 检查 fieldnames
        fnames = fieldnames(Models.PhasePipelineResult)
        @test :model_kind in fnames
        @test :xi in fnames
        @test :cep in fnames
        @test :first_order_boundary in fnames
        @test :spinodal in fnames
        @test :crossover_line in fnames
        @test :diagnostics in fnames
    end

    @testset "Production 类型存在" begin
        @test isdefined(Models, :FirstOrderSweepResult)
        @test isdefined(Models, :ProductionPipelineConfig)
        sweep_fields = fieldnames(Models.FirstOrderSweepResult)
        cfg_fields = fieldnames(Models.ProductionPipelineConfig)
        @test :first_point_fallback in sweep_fields
        @test :unknown_budget in cfg_fields
        @test :dT_initial in cfg_fields
    end

    @testset "run_phase_pipeline 接受临界区 direct re-evaluate 开关" begin
        tmp = mktempdir()
        result = Models.run_phase_pipeline(
            :PNJL;
            mode=:research,
            T_grid=[150.0],
            rho_grid=[0.1, 0.2, 0.3],
            xi=0.0,
            output_dir=tmp,
            profile=:smoke,
            solver_backend=:legacy,
            p_num=12,
            t_num=4,
            iterations=10,
            cep_strategy=:interpolate,
            cep_interpolate_use_direct_eval=true,
            promote_reference=false,
        )

        @test haskey(result.config_snapshot, "cep_interpolate_use_direct_eval")
        @test result.config_snapshot["cep_interpolate_use_direct_eval"] == true
    end

    @testset "production run_phase_pipeline 保留高精度 CEP 容差" begin
        tmp = mktempdir()
        result = Models.run_phase_pipeline(
            :PNJL;
            mode=:production,
            T_grid=[130.0, 131.0, 132.0],
            rho_grid=[0.0, 0.1, 0.2],
            xi=0.0,
            output_dir=tmp,
            profile=:smoke,
            solver_backend=:legacy,
            p_num=12,
            t_num=4,
            iterations=10,
            cep_tol=0.01,
            promote_reference=false,
        )

        @test haskey(result.config_snapshot, "cep_tol_MeV")
        @test result.config_snapshot["cep_tol_MeV"] == 0.01
    end

    @testset "run_phase_pipeline 暴露 research 与 production 模式" begin
        tmp_research = mktempdir()
        result_research = Models.run_phase_pipeline(
            :PNJL;
            mode=:research,
            T_grid=[150.0],
            rho_grid=[0.1, 0.2, 0.3],
            xi=0.0,
            output_dir=tmp_research,
            profile=:smoke,
            solver_backend=:legacy,
            p_num=12,
            t_num=4,
            iterations=10,
            promote_reference=false,
        )

        tmp_production = mktempdir()
        result_production = Models.run_phase_pipeline(
            :PNJL;
            mode=:production,
            T_grid=[150.0],
            rho_grid=[0.1, 0.2, 0.3],
            xi=0.0,
            output_dir=tmp_production,
            profile=:smoke,
            solver_backend=:legacy,
            p_num=12,
            t_num=4,
            iterations=10,
            promote_reference=false,
        )

        @test result_research.config_snapshot["mode"] == "research"
        @test result_production.config_snapshot["mode"] == "production"
    end

    @testset "run_phase_pipeline 非法 mode 抛出 ArgumentError" begin
        tmp = mktempdir()
        @test_throws ArgumentError Models.run_phase_pipeline(
            :PNJL;
            mode=:invalid,
            T_grid=[150.0],
            rho_grid=[0.1, 0.2, 0.3],
            xi=0.0,
            output_dir=tmp,
            profile=:smoke,
            solver_backend=:legacy,
            p_num=12,
            t_num=4,
            iterations=10,
            promote_reference=false,
        )
    end
end
