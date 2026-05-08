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
end
