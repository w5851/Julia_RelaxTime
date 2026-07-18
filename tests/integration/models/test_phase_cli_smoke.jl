using Test
using JSON3

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CLI_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "pnjl", "calculate_phase_structure.jl")

@testset "Phase CLI smoke" begin
    @test isfile(CLI_SCRIPT)

    output_dir = mktempdir()
    cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) --preset=smoke --iterations=6 --p_num=8 --t_num=4 --cep_max_bisect_iter=1 --cep_max_refine_level=0 --output_dir=$(output_dir)`
    run(cmd)

    summary_path = joinpath(output_dir, "phase_summary.json")
    manifest_path = joinpath(output_dir, "run_manifest.json")
    grid_convergence_path = joinpath(output_dir, "phase_grid_convergence.csv")

    @test isfile(summary_path)
    @test isfile(manifest_path)
    @test isfile(grid_convergence_path)

    summary = JSON3.read(read(summary_path, String))
    manifest = JSON3.read(read(manifest_path, String))

    @test haskey(summary, "stats")
    @test haskey(manifest, "preset")
    @test String(manifest["preset"]) == "smoke"
    @test haskey(manifest, "effective_config")
    @test haskey(manifest["effective_config"], "rho_geometry_convergence")
end
