using Test

const _MUB_SENS_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(_MUB_SENS_ROOT, "scripts", "analysis", "relaxtime", "scan_charged_gbu_mub_sensitivity.jl"))
const _MUB_SENS = ChargedGBUMuBSensitivity

@testset "charged GBU muB sensitivity algebra and contract" begin
    @test _MUB_SENS.parse_float_list("3,5,7.7", "energies") == [3.0, 5.0, 7.7]
    @test_throws ArgumentError _MUB_SENS.parse_float_list("3,3", "energies")
    @test_throws ArgumentError _MUB_SENS.validate_grids([3.0], [-50.0, 50.0])
    opts = _MUB_SENS.parse_args(["--sqrts-list", "3,7.7", "--mub-offsets", "-50,0,50",
        "--q-nodes", "8", "--mesh", "32", "--no-plot"])
    @test opts[:energies] == [3.0, 7.7] && opts[:offsets] == [-50.0, 0.0, 50.0]
    @test opts[:settings].q_nodes == 8 && opts[:settings].mesh == 32 && !opts[:make_plot]

    points = _MUB_SENS.probe_points([7.7], [-50.0, 0.0, 50.0])
    @test length(points) == 3
    @test [p.delta_muB_MeV for p in points] == [0.0, -50.0, 50.0]
    @test all(p.T_MeV == points[1].T_MeV for p in points)
    @test points[1].on_baseline_freezeout
    @test !points[2].on_baseline_freezeout && !points[3].on_baseline_freezeout

    channels = Dict{String,Any}(
        "pi_plus" => (passed=true, density_inv_fm3=2.0),
        "K_plus" => (passed=true, density_inv_fm3=1.0),
    )
    @test _MUB_SENS.ratio_from_channels(channels).ratio == 0.5
    channels["pi_plus"] = (passed=false, density_inv_fm3=2.0)
    @test !_MUB_SENS.ratio_from_channels(channels).passed

    rows = [
        (sqrt_s_NN_GeV=3.0, delta_muB_MeV=-50.0, ratio=0.40),
        (sqrt_s_NN_GeV=3.0, delta_muB_MeV=0.0, ratio=0.50),
        (sqrt_s_NN_GeV=3.0, delta_muB_MeV=50.0, ratio=0.60),
    ]
    d = _MUB_SENS.sensitivity_rows(rows; delta_muB_MeV=50.0)
    @test length(d) == 1 && d[1].passed
    @test d[1].d_ratio_dmuB ≈ 0.002
    @test d[1].d_logratio_dmuB ≈ (log(0.60) - log(0.40)) / 100.0
    @test d[1].curvature ≈ 0.0

    @test occursin("BQS re-solved", read(joinpath(_MUB_SENS_ROOT,
        "scripts", "analysis", "relaxtime", "scan_charged_gbu_mub_sensitivity.jl"), String))
    @test occursin("production_authorized=false", read(joinpath(_MUB_SENS_ROOT,
        "scripts", "analysis", "relaxtime", "scan_charged_gbu_mub_sensitivity.jl"), String))
end
