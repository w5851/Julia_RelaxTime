# TmuScan.jl 单元测试
#
# 测试内容：
# 1. 模块加载与常量
# 2. run_tmu_scan 接口存在

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

const _TMS = Models.TmuScan

# ============================================================================

@testset "TmuScan" begin

    @testset "DEFAULT_T_VALUES 非空" begin
        @test length(_TMS.DEFAULT_T_VALUES) > 0
        @test all(t -> t > 0 && isfinite(t), _TMS.DEFAULT_T_VALUES)
    end

    @testset "DEFAULT_MU_VALUES 非空" begin
        @test length(_TMS.DEFAULT_MU_VALUES) > 0
        @test all(isfinite, _TMS.DEFAULT_MU_VALUES)
    end

    @testset "DEFAULT_OUTPUT_PATH 非空" begin
        @test !isempty(_TMS.DEFAULT_OUTPUT_PATH)
    end

    @testset "run_tmu_scan 接口存在" begin
        @test isdefined(_TMS, :run_tmu_scan)
        @test _TMS.run_tmu_scan isa Function
    end

    @testset "auto backend 语义路由配置" begin
        @test _TMS._effective_solver_backend(:auto, :PNJL; auto_pnjl_backend=:models) == :models
        @test _TMS._effective_solver_backend(:auto, :PNJL; auto_pnjl_backend=:legacy) == :legacy
        @test _TMS._effective_solver_backend(:auto, :RPNJL; auto_pnjl_backend=:legacy) == :models
    end

    @testset "semantic_mode 参数校验" begin
        @test _TMS._validate_semantic_mode(:ground_state, nothing) === nothing
        @test_throws ArgumentError _TMS._validate_semantic_mode(:constrained_manifold, nothing)
        @test_throws ArgumentError _TMS._validate_semantic_mode(:ground_state, (_h, _q) -> :hadron)
        @test_throws ArgumentError _TMS._validate_semantic_mode(:invalid_mode, nothing)
    end

    @testset "models 路径禁止 legacy solver 开关" begin
        @test _TMS._reject_legacy_solver_kwargs((; solver=:newton)) === nothing
        @test_throws ArgumentError _TMS._reject_legacy_solver_kwargs((; use_problem_spec=false))
        @test_throws ArgumentError _TMS._reject_legacy_solver_kwargs((; allow_legacy_path=true))
        @test_throws ArgumentError _TMS._reject_legacy_solver_kwargs((; warn_on_legacy_path=false))
        @test_throws ArgumentError _TMS._reject_legacy_solver_kwargs((; problem_spec=:dummy))
    end
end
