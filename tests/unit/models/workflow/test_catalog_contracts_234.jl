using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const EXPECTED_STAGE_IDS = (
    :prepare_inputs,
    :solve_core,
    :postprocess,
    :export_artifacts,
    :emit_diagnostics,
    :emit_repro_manifest,
)

@testset "catalog contracts 2/3/4" begin
    @testset "workflow/scan/orchestrator share ordered stage skeleton" begin
        @test Models.workflow_pipeline_stage_ids() == EXPECTED_STAGE_IDS
        @test Models.scan_pipeline_stage_ids() == EXPECTED_STAGE_IDS
        @test Models.relaxtime_orchestrator_stage_ids() == EXPECTED_STAGE_IDS
    end
end
