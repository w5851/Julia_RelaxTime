using Test
using JSON3

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "run_phase_guided_transport_scan.jl")

@testset "phase guided transport scan dry-run smoke" begin
    @test isfile(SCRIPT_PATH)

    outdir = mktempdir()
    cmd = `julia --project=. $SCRIPT_PATH --mode fixed-muB-phase-scaled --outdir $outdir --case-name smoke_case --xi-list -0.2,0.0,0.2 --muB-list 150,400 --alphaT-list 1.0,1.1 --dry-run --overwrite`
    run(cmd)

    plan_csv = joinpath(outdir, "sampling_plan.csv")
    readme = joinpath(outdir, "README.md")
    eff_cfg = joinpath(outdir, "effective_config.json")
    manifest = joinpath(outdir, "run_manifest.json")

    @test isfile(plan_csv)
    @test isfile(readme)
    @test isfile(eff_cfg)
    @test isfile(manifest)

    plan_text = read(plan_csv, String)
    @test occursin("phase_reference_kind", plan_text)
    @test occursin("mode_a_fixed_muB_phase_scaled", plan_text)

    readme_text = read(readme, String)
    @test occursin("phase-guided transport canonical case", readme_text)
    @test occursin("sampling_plan.csv", readme_text)
    @test occursin("smoke_case", readme_text)

    cfg_text = read(eff_cfg, String)
    @test occursin("\"mode\":\"mode_a_fixed_muB_phase_scaled\"", cfg_text)

    manifest_obj = JSON3.read(read(manifest, String))
    @test String(manifest_obj["script"]) == "scripts/relaxtime/run_phase_guided_transport_scan.jl"
end
