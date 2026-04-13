using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const ORCH_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "run_relaxtime_orchestrator.jl")

function _write_lightweight_cfg(path::String)
    open(path, "w") do io
        write(io, """
schema_version = "v1"

[scan.transport]
muB_MeV = 0.0
xi_list = [0.0]
tmin_MeV = 150.0
tmax_MeV = 150.0
tstep_MeV = 10.0
resume = false
overwrite = true

[scan.transport.solver]
p_num = 8
t_num = 2
max_iter = 20

[scan.transport.tau]
mode = "finite_15"
tau_p_nodes = 8
tau_angle_nodes = 2
tau_phi_nodes = 2
tau_n_sigma = 4
sigma_grid_n = 32

[scan.transport.transport]
compute_bulk = false
tr_p_nodes = 12
tr_p_max_fm = 8.0

[scan.cross_section]
muB_MeV = 0.0
T_list_MeV = [150.0]
xi_list = [0.0]
processes = ["ud_to_ud"]
n_points = 16

[scan.cross_section.energy]
mode = "linspace"
sqrt_s_min_MeV = 400.0
sqrt_s_max_MeV = 1000.0
sqrt_s_num = 16
""")
    end
end

@testset "relaxtime orchestrator smoke" begin
    @test isfile(ORCH_PATH)

    outdir = mktempdir()
    cfg = joinpath(outdir, "orchestrator_smoke_cfg.toml")
    _write_lightweight_cfg(cfg)

    cmd_transport = `julia --project=. $ORCH_PATH transport --config $cfg --outdir $outdir --resume`
    run(cmd_transport)
    @test isfile(joinpath(outdir, "run_manifest.json"))
    @test isfile(joinpath(outdir, "effective_config.json"))
    @test isfile(joinpath(outdir, "consumption_report.json"))

    cmd_xs = `julia --project=. $ORCH_PATH cross-section --config $cfg --outdir $outdir --processes udtoud`
    run(cmd_xs)
    @test isfile(joinpath(outdir, "run_manifest.json"))
    @test isfile(joinpath(outdir, "effective_config.json"))
    @test isfile(joinpath(outdir, "consumption_report.json"))

    fallback_file = joinpath(outdir, "fallback_events.json")
    open(fallback_file, "w") do io
        write(io, "[{\"event\":\"bulk_fallback\",\"point\":\"T=150,muB=0,xi=0.0\"}]")
    end
    cmd_fail = `julia --project=. $ORCH_PATH transport --config $cfg --outdir $outdir --fail-on-fallback`
    @test_throws Base.ProcessFailedException run(cmd_fail)
end
