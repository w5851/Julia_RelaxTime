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

    @testset "CSV field escaping" begin
        @test Models._phase_csv_value(nothing) == ""
        @test Models._phase_csv_value("plain") == "plain"
        @test Models._phase_csv_value("a,b") == "\"a,b\""
        @test Models._phase_csv_value("a\"b") == "\"a\"\"b\""
        @test Models._phase_csv_value("a\nb") == "\"a\nb\""
        @test Models._phase_csv_value("a\rb") == "\"a\rb\""
    end

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
        @test haskey(paths, "phase_grid_convergence")
        @test isfile(paths["phase_grid_convergence"])
        @test startswith(read(paths["phase_grid_convergence"], String), "axis,xi,T_MeV")
    end

    @testset "grid convergence reasons are valid quoted CSV fields" begin
        tmp = mktempdir()
        reason = "midpoint_classification_changed_or_unresolved:valid,unknown,\"valid\"\nreview"
        result = Models.PhasePipelineResult(
            xi=0.3,
            diagnostics=Dict{String, Any}(
                "grid_convergence_records" => [(
                    axis="temperature",
                    xi=0.3,
                    T_MeV=10.0,
                    level=1,
                    left=5.0,
                    right=15.0,
                    midpoint=10.0,
                    position_error_MeV=0.1,
                    density_error=0.01,
                    maxwell_area=1e-4,
                    response_rtol=0.05,
                    converged=false,
                    reason=reason,
                )],
            ),
        )
        paths = Models.build_phase_artifacts(result; output_dir=tmp)
        grid_csv = read(paths["phase_grid_convergence"], String)
        @test occursin(
            "\"midpoint_classification_changed_or_unresolved:valid,unknown,\"\"valid\"\"\nreview\"",
            grid_csv,
        )
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
