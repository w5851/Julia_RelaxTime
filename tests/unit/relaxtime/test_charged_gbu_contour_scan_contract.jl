using Test, JSON3

const CONTOUR_SCAN_PATH = joinpath(@__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime", "run_charged_gbu_contour_scan.jl")
isdefined(Main, :ChargedGBUContourScan) || include(CONTOUR_SCAN_PATH)
const CGBU_CONTOUR_SCAN = Main.ChargedGBUContourScan

@testset "Charged GBU contour scanner pure contracts" begin
    source = read(CONTOUR_SCAN_PATH, String)
    @test occursin("charged_gbu_infinite_qgate", source) || occursin("Workflow.I.kernel", source)
    @test occursin("source_hashes", source)
    @test occursin("background=bg === nothing ? nothing : (residual=bg.residual, elapsed_s=elapsed_s, seed=seed)", source)
    @test !occursin("CausalGBUResearch", source)
    @test CGBU_CONTOUR_SCAN.parse_grid("40:220:10") == collect(40.0:10.0:220.0)
    @test CGBU_CONTOUR_SCAN.parse_grid("0:800:50")[end] == 800.0
    @test_throws ArgumentError CGBU_CONTOUR_SCAN.parse_grid("40:220:7")
    @test_throws ArgumentError CGBU_CONTOUR_SCAN.parse_grid("40:220:0")
    @test CGBU_CONTOUR_SCAN.parse_channels("pi_plus,K_plus") == [:pi_plus, :K_plus]
    @test_throws ArgumentError CGBU_CONTOUR_SCAN.parse_channels("pi_plus,pi_plus")
    @test_throws ArgumentError CGBU_CONTOUR_SCAN.parse_channels("pi_zero")
    @test CGBU_CONTOUR_SCAN.shard_indices(7, 0, 3) == [1, 4, 7]
    @test CGBU_CONTOUR_SCAN.shard_indices(7, 2, 3) == [3, 6]
    @test CGBU_CONTOUR_SCAN.point(80, 719).muq_MeV == 719 / 3
    opts = CGBU_CONTOUR_SCAN.parse_args(["--output", "x", "--t-grid", "40:60:10", "--resume"])
    @test opts[:resume] && opts[:output] == "x"
    @test CGBU_CONTOUR_SCAN.parse_args(["--help"]) === nothing
    @test opts[:density_route] == "direct_finite_q"
    @test CGBU_CONTOUR_SCAN.parse_args(["--output", "x", "--density-route", "q0_lambda_reference"])[:density_route] == "q0_lambda_reference"
    @test_throws ArgumentError CGBU_CONTOUR_SCAN.parse_args(["--output", "x", "--density-route", "folded"])
    settings = CGBU_CONTOUR_SCAN.screening_settings()
    @test settings == (mesh=64, cut_nodes=32, tail_nodes=32, omega_nodes=64, q_nodes=8, qmax=8.)
    direct = CGBU_CONTOUR_SCAN._scan_identity(Dict(), [40.], [0.], [:pi_plus], settings, 0, 1, Dict())
    extrapolated = CGBU_CONTOUR_SCAN._scan_identity(Dict(), [40.], [0.], [:pi_plus], settings, 0, 1, Dict(); route="q0_lambda_reference")
    @test direct != extrapolated # Never resume a different algorithm's points.
    @test occursin("external omega", CGBU_CONTOUR_SCAN.coordinate_contract("q0_lambda_reference").bose_frequency)
    finite = CGBU_CONTOUR_SCAN.screening_settings(route="q0_lambda_reference")
    @test finite.endpoint_policy == "finite_window" && finite.omega_lower_inv_fm == 1e-5
    fine = merge(finite, (omega_lower_inv_fm=1e-6,))
    id(s; bg=nothing) = CGBU_CONTOUR_SCAN._scan_identity(Dict(), [40.], [0.], [:pi_plus], s, 0, 1, Dict();
        route="q0_lambda_reference", background_fingerprint=bg)
    @test id(finite) != id(fine)
    @test id(finite) != id(merge(finite, (endpoint_policy="strict_zero_limit",)))
    @test id(finite; bg="savedA") != id(finite; bg="savedB")
    opts = CGBU_CONTOUR_SCAN.parse_args(["--output", "x", "--q0-omega-lower", "1e-6", "--background-input-root", "saved"])
    @test opts[:q0_omega_lower] == 1e-6 && opts[:background_input_root] == "saved"
    @test_throws ArgumentError CGBU_CONTOUR_SCAN.parse_args(["--output", "x", "--q0-omega-lower", "0"])
end

@testset "Frozen contour backgrounds: source, grid and resume identity" begin
    s = CGBU_CONTOUR_SCAN; b = s.SavedBackgrounds
    model = Main.Models.create_model(:PNJL)
    seed = [-1., -1.1, -1.4, .1, .2, .5, .6, .2]
    restored = b.restore_seed(model, seed, 140., 425., 1e-14)
    @test restored.mu == (u=.5, d=.6, s=.2)
    @test restored.Phi == .1 && restored.PhiBar == .2
    @test collect(restored.m) == Main.Models.calculate_mass_vec(model, Main.Models.meanfield_state(seed[1:5]).phi)
    @test_throws ArgumentError b.restore_seed(model, seed[1:7], 140., 425., 1e-14)
    mktempdir() do dir
        manifest = (schema="charged_gbu_contour_scan_v2", T_grid=[140.], muB_grid=[425.],
            scan_identity="fixture", git_head="fixture", source_hashes=s.source_hashes(s.DEFAULT_CONFIG))
        s._write_atomic(joinpath(dir,"manifest.json"), manifest)
        point = (schema="charged_gbu_contour_point_v2", scan_identity="fixture", row_index=1,
            col_index=1, T_MeV=140., muB_MeV=425., background=(seed=seed, residual=1e-14))
        path = joinpath(dir,"points","point.json"); s._write_atomic(path,point)
        first = b.read_snapshot(dir, [140.], [425.])
        @test !first.provenance.solver_called && length(first.points) == 1
        @test_throws ArgumentError b.read_snapshot(dir, [145.], [425.])
        s._write_atomic(path,merge(point,(background=(seed=2seed,residual=1e-14),)))
        @test first.provenance.fingerprint != b.read_snapshot(dir,[140.],[425.]).provenance.fingerprint
        s._write_atomic(path,merge(point,(scan_identity="wrong",)))
        @test_throws ArgumentError b.read_snapshot(dir,[140.],[425.])
    end
end
