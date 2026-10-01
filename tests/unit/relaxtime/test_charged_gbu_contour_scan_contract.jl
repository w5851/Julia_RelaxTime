using Test

const CONTOUR_SCAN_PATH = joinpath(@__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime", "run_charged_gbu_contour_scan.jl")
include(CONTOUR_SCAN_PATH)
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
end
