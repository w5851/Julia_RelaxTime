# PhaseArtifacts.jl 单元测试
#
# 测试内容：
# 1. resolve_phase_output_target 路径构建
# 2. build_phase_artifacts / promote_phase_artifacts 接口存在
# 3. PhasePipelineResult 结构

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

# ============================================================================

@testset "PhaseArtifacts" begin

    @testset "resolve_phase_output_target 接口" begin
        @test isdefined(Models, :resolve_phase_output_target)
        target = Models.resolve_phase_output_target(:PNJL; project_root=PROJECT_ROOT)
        @test target isa NamedTuple
    end

    @testset "build_phase_artifacts 接口存在" begin
        @test isdefined(Models, :build_phase_artifacts)
        @test Models.build_phase_artifacts isa Function
    end

    @testset "promote_phase_artifacts 接口存在" begin
        @test isdefined(Models, :promote_phase_artifacts)
        @test Models.promote_phase_artifacts isa Function
    end

    @testset "PhasePipelineResult 结构" begin
        @test isdefined(Models, :PhasePipelineResult)
    end

    @testset "PromotionResult 结构" begin
        @test isdefined(Models, :PromotionResult)
        res = Models.PromotionResult(
            passed=true,
            baseline_id="test",
            failed_checks=String[],
            reference_dir="",
        )
        @test res.passed == true
        @test isempty(res.failed_checks)
    end
end
