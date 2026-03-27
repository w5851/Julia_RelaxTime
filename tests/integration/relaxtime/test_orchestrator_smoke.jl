using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const ORCH_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "run_relaxtime_orchestrator.jl")

@testset "relaxtime orchestrator smoke" begin
    @test isfile(ORCH_PATH)

    outdir = mktempdir()
    cfg = joinpath(REPO_ROOT, "config", "workflows", "relaxtime", "default.toml")

    cmd_transport = `julia --project=. $ORCH_PATH transport --config $cfg --outdir $outdir --resume`
    run(cmd_transport)
    @test isfile(joinpath(outdir, "run_manifest.json"))
    @test isfile(joinpath(outdir, "effective_config.json"))
    @test isfile(joinpath(outdir, "consumption_report.json"))

    cmd_xs = `julia --project=. $ORCH_PATH cross-section --config $cfg --outdir $outdir --processes udtoud,ustous`
    run(cmd_xs)
    @test isfile(joinpath(outdir, "run_manifest.json"))
    @test isfile(joinpath(outdir, "effective_config.json"))
    @test isfile(joinpath(outdir, "consumption_report.json"))
end
