using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const ORCH_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "run_relaxtime_orchestrator.jl")

@testset "plot contract smoke" begin
    @test isfile(ORCH_PATH)

    outdir = mktempdir()
    cfg = joinpath(outdir, "plot_cfg.toml")
    open(cfg, "w") do io
        write(io, """
schema_version = "v1"

[scan.transport]
muB_MeV = 0.0

[scan.cross_section]
muB_MeV = 0.0
T_list_MeV = [150.0, 250.0]
processes = ["ud_to_ud", "us_to_us"]

[plot.transport]
ys = ["eta_over_s", "sigma_over_T"]
""")
    end

    run(`julia --project=. $ORCH_PATH transport --config $cfg --outdir $outdir`)
    run(`julia --project=. $ORCH_PATH cross-section --config $cfg --outdir $outdir`)

    fig_dir = joinpath(outdir, "figures")
    @test isdir(fig_dir)

    expected_transport = Set([
        "transport__eta_over_s__muB0.png",
        "transport__sigma_over_T__muB0.png",
    ])
    for f in expected_transport
        @test isfile(joinpath(fig_dir, f))
    end

    expected_xs = Set([
        "xsec__T150__ud_to_ud.png",
        "xsec__T150__us_to_us.png",
        "xsec__T250__ud_to_ud.png",
        "xsec__T250__us_to_us.png",
    ])
    for f in expected_xs
        @test isfile(joinpath(fig_dir, f))
    end
end
