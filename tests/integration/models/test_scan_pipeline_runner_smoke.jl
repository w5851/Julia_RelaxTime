using Test
using JSON3

_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

@testset "Scan pipeline runner smoke" begin
    tmp = mktempdir()
    output_path = joinpath(tmp, "scan_tmu.csv")

    result = Models.run_scan_pipeline(
        :tmu;
        model_kind=:PNJL,
        T_values=[150.0],
        mu_values=[0.0],
        xi_values=[0.0],
        output_path=output_path,
        overwrite=true,
        resume=false,
        use_phase_aware=false,
        solver_backend=:models,
        p_num=12,
        t_num=4,
        iterations=40,
    )

    @test hasproperty(result, :manifest_path)
    @test isfile(result.manifest_path)
    @test isfile(output_path)

    manifest = JSON3.read(read(result.manifest_path, String))
    @test haskey(manifest, :pipeline)
    @test String(manifest.pipeline.pipeline_family) == "scan"
end
