using Test
using JSON3

_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

@testset "Workflow adapter pipeline smoke" begin
    tmp = mktempdir()

    result = Models.run_workflow_pipeline(
        :transport;
        T_fm=0.15,
        mu_fm=0.0,
        xi=0.0,
        output_dir=tmp,
        compute_tau=false,
        tau=(u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0),
        compute_bulk=false,
    )

    @test hasproperty(result, :run_manifest)
    manifest_path = result.run_manifest
    @test isfile(manifest_path)

    manifest = JSON3.read(read(manifest_path, String))
    @test haskey(manifest, :pipeline)
    @test String(manifest.pipeline.pipeline_family) == "workflow"
end
