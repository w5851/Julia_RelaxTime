using Test
using JSON3

_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

const EXPECTED_STAGE_SEQUENCE = [
    "build_model",
    "prepare_grid",
    "solve_points",
    "collect_diagnostics",
    "analyze_phase",
    "export_artifacts",
    "emit_repro_manifest",
]

@testset "Phase pipeline runner smoke" begin
    tmp = mktempdir()

    result = Models.run_phase_pipeline(
        :PNJL;
        mode=:research,
        T_grid=[150.0],
        rho_grid=[0.1, 0.2, 0.3],
        xi=0.0,
        output_dir=tmp,
        profile=:smoke,
        solver_backend=:models,
        p_num=12,
        t_num=4,
        iterations=80,
        promote_reference=false,
    )

    @test haskey(result.artifact_paths, "phase_summary")
    @test haskey(result.artifact_paths, "phase_report")
    @test haskey(result.diagnostics, "scan_total")

    manifest_path = joinpath(tmp, "run_manifest.json")
    @test isfile(manifest_path)

    manifest = JSON3.read(read(manifest_path, String))
    @test haskey(manifest, :completed_stages)
    @test collect(String.(manifest.completed_stages)) == EXPECTED_STAGE_SEQUENCE
    @test haskey(manifest, :pipeline)
    @test haskey(manifest.pipeline, :artifact_hash)
    @test !isempty(String(manifest.pipeline.artifact_hash))
end
