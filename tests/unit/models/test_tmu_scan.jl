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
end
