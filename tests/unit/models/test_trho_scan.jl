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
end
