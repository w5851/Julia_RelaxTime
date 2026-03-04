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
end
