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

    @testset "CEP artifact exposes explicit mu contracts" begin
        tmp = mktempdir()
        result = Models.PhasePipelineResult(
            cep=Models.CEPResult(found=true, T_cep_MeV=130.0, mu_cep_MeV=292.0),
        )
        paths = Models.build_phase_artifacts(result; output_dir=tmp)
        summary = read(paths["phase_summary"], String)
        report = read(paths["phase_report"], String)

        @test occursin("\"muq_cep_MeV\":292", summary)
        @test occursin("\"muB_cep_MeV\":876", summary)
        @test occursin("\"mu_cep_MeV\":292", summary)
        @test occursin("- muq_cep_MeV: 292.0", report)
        @test occursin("- muB_cep_MeV: 876.0", report)
        @test occursin("compatibility alias for muq_cep_MeV", report)
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
