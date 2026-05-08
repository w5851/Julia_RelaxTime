using Test
using JSON3

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCAN_SCRIPT = joinpath(REPO_ROOT, "scripts", "relaxtime", "run_phase_guided_transport_scan.jl")
const PLOT_SCRIPT = joinpath(REPO_ROOT, "scripts", "relaxtime", "run_phase_guided_transport_plots.jl")

@testset "phase guided transport plot smoke" begin
    python_cmd = Sys.which("python") !== nothing ? "python" : (Sys.which("python3") !== nothing ? "python3" : nothing)
    has_matplotlib = python_cmd !== nothing && success(`$(python_cmd) -c "import matplotlib"`)
    if !has_matplotlib
        @test true
        return
    end

    outdir = mktempdir()
    fig_dir = mktempdir()
    run(`julia --project=. $SCAN_SCRIPT --mode fixed-T-sparse-muB --outdir $outdir --case-name plot_smoke --xi-list -0.2,0.0,0.2 --muB-list 0,450 --T-list 120 --overwrite`)
    run(`julia --project=. $PLOT_SCRIPT --case-dir $outdir --fig-dir $fig_dir --python $python_cmd --overwrite`)

    manifest_path = joinpath(fig_dir, "plot_manifest.json")
    readme_path = joinpath(outdir, "README.md")

    @test isdir(fig_dir)
    @test isfile(manifest_path)
    @test isfile(readme_path)

    manifest = JSON3.read(read(manifest_path, String))
    @test Int(manifest["count"]) >= 8

    readme = read(readme_path, String)
    @test occursin("Generated Figures", readme)
    @test occursin("plot_manifest.json", readme)
    @test occursin(replace(fig_dir, '\\' => '/'), readme)
    @test occursin("zeta_vs_xi.png", readme)
end
