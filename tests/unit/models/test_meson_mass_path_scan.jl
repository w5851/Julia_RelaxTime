using Test

const PROJECT_ROOT_MMPS = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_MMPS, "src", "models", "Models.jl"))
end

const _MMPS = Models.MesonMassPathScan

@testset "MesonMassPathScan" begin
    @testset "default output paths exist as strings" begin
        @test _MMPS.DEFAULT_MESON_MASS_PATH_OUTPUT_PATH isa String
        @test _MMPS.DEFAULT_FREEZEOUT_MESON_MASS_OUTPUT_PATH isa String
        @test _MMPS.DEFAULT_ISENTROPIC_MESON_MASS_OUTPUT_PATH isa String
    end

    @testset "entrypoints exist" begin
        @test isdefined(_MMPS, :run_meson_mass_path_scan)
        @test isdefined(_MMPS, :run_freezeout_meson_mass_scan)
        @test isdefined(_MMPS, :run_isentropic_meson_mass_scan)
        @test Models.run_meson_mass_path_scan isa Function
        @test Models.run_freezeout_meson_mass_scan isa Function
        @test Models.run_isentropic_meson_mass_scan isa Function
    end

    @testset "header keeps path metadata first" begin
        @test startswith(_MMPS.HEADER, "path_family,path_profile,path_segment,path_point_index,path_order_key,path_label,")
        @test occursin("equilibrium_residual_norm", _MMPS.HEADER)
        @test occursin("M_pi", _MMPS.HEADER)
        @test occursin("Gamma_sigma_prime", _MMPS.HEADER)
    end

    @testset "point normalization distinguishes fixed-mu and fixed-sigma paths" begin
        fixedmu = _MMPS._normalize_path_point((T_MeV=160.0, muB_MeV=90.0, path_family="freezeout"), 1)
        fixedsigma = _MMPS._normalize_path_point((T_MeV=160.0, sigma_target=25.0), 1)
        @test fixedmu.path_family == "freezeout"
        @test isfinite(fixedmu.muB_MeV)
        @test fixedsigma.path_family == "isentropic"
        @test !isfinite(fixedsigma.muB_MeV)
        @test fixedsigma.sigma_target == 25.0
    end

    @testset "row formatting keeps integer index and equilibrium iterations" begin
        cols = ["path_point_index", "equilibrium_iterations", "message"]
        row = Dict{String,Any}(
            "path_point_index" => 3,
            "equilibrium_iterations" => -1,
            "message" => "",
        )
        values = _MMPS._row_to_values(cols, row)
        @test values[1] == "3"
        @test values[2] == "-1"
    end

    @testset "schedule points preserves input order" begin
        points = [
            (path_family="freezeout", path_profile="p", path_label="a", path_order_key=3.0, path_point_index=0),
            (path_family="freezeout", path_profile="p", path_label="a", path_order_key=1.0, path_point_index=1),
        ]
        scheduled = _MMPS._schedule_points(copy(points))
        @test scheduled == points
    end
end
