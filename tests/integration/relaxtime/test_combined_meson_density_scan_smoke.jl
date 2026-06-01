using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT = joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_combined_meson_density_scan.jl")
const OUTDIR = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "scan", "test_outputs", "combined_meson_density")
const FIGDIR = joinpath(PROJECT_ROOT, "data", "outputs", "figures", "relaxtime", "scan", "test_outputs", "combined_meson_density")
const OUTCSV = joinpath(OUTDIR, "combined_meson_density_scan.csv")
const OUTSVG = joinpath(FIGDIR, "combined_meson_density_scan.svg")
const PLOTMANIFEST = joinpath(FIGDIR, "plot_manifest.json")
const OUTREADME = joinpath(OUTDIR, "README.md")

@testset "combined meson density scan smoke" begin
    isdir(OUTDIR) && rm(OUTDIR; recursive=true, force=true)
    isdir(FIGDIR) && rm(FIGDIR; recursive=true, force=true)
    mkpath(OUTDIR)

    cmd = `julia --project=. $SCRIPT --output-dir $OUTDIR --overwrite --tmin 150 --tmax 150 --tstep 10 --muq 0 --p-num 4 --t-num 2 --max-iter 6 --stable-q-nodes 6 --qmax 2 --q-nodes 3 --omega-max 2 --omega-nodes 3 --phase-display fold_0_pi`
    run(cmd)

    @test isfile(OUTCSV)
    @test isfile(OUTSVG)
    @test isfile(PLOTMANIFEST)
    @test isfile(OUTREADME)

    text = read(OUTCSV, String)
    @test occursin("# bridge: path_strategy x density_regime", text)
    @test occursin("path_strategy,path_point_index,T_MeV", text)
    @test occursin("phase_display", text)
    @test occursin("stable", text)
    @test occursin("strict_bw_stage1", text)
    @test occursin("phase_shift_current", text)
    @test occursin("phase_shift_gbu_reference", text)

    data_lines = [
        line for line in split(text, '\n')
        if !isempty(strip(line)) && !startswith(strip(line), "#") && !startswith(strip(line), "path_strategy,")
    ]
    @test length(data_lines) == 4

    summary = read(OUTREADME, String)
    @test occursin("Bridge-style composition", summary)
    @test occursin("figure directory", summary)
    @test occursin("Status Counts", summary)

    svg = read(OUTSVG, String)
    @test occursin("<svg", svg)
    @test occursin("Combined Meson Density Scan", svg)

    manifest = read(PLOTMANIFEST, String)
    @test occursin("combined_meson_density_plot_manifest_v1", manifest)
    @test occursin("combined_meson_density_scan.svg", manifest)
end
