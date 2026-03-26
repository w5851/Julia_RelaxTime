using Test
using JSON3

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CLI_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "pnjl", "calculate_phase_structure.jl")

@testset "Phase CLI smoke" begin
    @test isfile(CLI_SCRIPT)

    output_dir = mktempdir()
    cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) --preset=smoke --output_dir=$(output_dir)`
    run(cmd)

    summary_path = joinpath(output_dir, "phase_summary.json")
    manifest_path = joinpath(output_dir, "run_manifest.json")

    @test isfile(summary_path)
    @test isfile(manifest_path)

    summary = JSON3.read(read(summary_path, String))
    manifest = JSON3.read(read(manifest_path, String))

    @test haskey(summary, "cep")
    @test haskey(summary, "stats")
    @test haskey(manifest, "preset")
    @test String(manifest["preset"]) == "smoke"
    @test haskey(manifest, "effective_config")
end
