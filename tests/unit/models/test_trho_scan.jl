# TrhoScan.jl 单元测试
#
# 测试内容：
# 1. 模块加载与常量
# 2. build_default_rho_grid
# 3. run_trho_scan 接口存在

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# TrhoScan 是 Models 的子模块
const _TS = Models.TrhoScan

# ============================================================================

@testset "TrhoScan" begin

    @testset "DEFAULT_T_VALUES 非空" begin
        @test length(_TS.DEFAULT_T_VALUES) > 0
        @test all(t -> t > 0 && isfinite(t), _TS.DEFAULT_T_VALUES)
    end

    @testset "DEFAULT_RHO_VALUES 非空" begin
        @test length(_TS.DEFAULT_RHO_VALUES) > 0
        @test all(r -> isfinite(r), _TS.DEFAULT_RHO_VALUES)
    end

    @testset "DEFAULT_OUTPUT_PATH 非空" begin
        @test !isempty(_TS.DEFAULT_OUTPUT_PATH)
    end

    @testset "build_default_rho_grid" begin
        grid = _TS.build_default_rho_grid()
        @test grid isa AbstractVector
        @test length(grid) > 0
        @test issorted(grid)
    end

    @testset "run_trho_scan 接口存在" begin
        @test isdefined(_TS, :run_trho_scan)
        @test _TS.run_trho_scan isa Function
    end

    @testset "auto backend 语义路由配置" begin
        @test _TS._effective_solver_backend(:auto, :PNJL; auto_pnjl_backend=:models) == :models
        @test _TS._effective_solver_backend(:auto, :PNJL; auto_pnjl_backend=:legacy) == :legacy  # route preserved; execution guard rejects legacy backend
        @test _TS._effective_solver_backend(:auto, :RPNJL; auto_pnjl_backend=:legacy) == :models
    end

    @testset "legacy backend removed from input validation" begin
        @test_throws ArgumentError _TS._validate_trho_scan_inputs([150.0], [0.2], [0.0], :candidates, :fixed_rho, :legacy, :PNJL)
    end

    @testset "semantic_mode 参数校验" begin
        @test _TS._validate_semantic_mode(:ground_state) === nothing
        @test _TS._validate_semantic_mode(:constrained_manifold) === nothing
        @test_throws ArgumentError _TS._validate_semantic_mode(:invalid_mode)
    end

    @testset "auto_pnjl_backend 参数校验" begin
        @test _TS._validate_auto_pnjl_backend(:models) === nothing
        @test _TS._validate_auto_pnjl_backend(:legacy) === nothing
        @test_throws ArgumentError _TS._validate_auto_pnjl_backend(:invalid_backend)
    end

    @testset "models 路径禁止 legacy solver 开关" begin
        @test _TS._reject_legacy_solver_kwargs((; solver=:newton)) === nothing
        @test_throws ArgumentError _TS._reject_legacy_solver_kwargs((; use_problem_spec=false))
        @test_throws ArgumentError _TS._reject_legacy_solver_kwargs((; allow_legacy_path=true))
        @test_throws ArgumentError _TS._reject_legacy_solver_kwargs((; warn_on_legacy_path=false))
        @test_throws ArgumentError _TS._reject_legacy_solver_kwargs((; problem_spec=:dummy))
    end
end
