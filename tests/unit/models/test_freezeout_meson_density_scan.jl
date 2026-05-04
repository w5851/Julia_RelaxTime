using Test

const PROJECT_ROOT_FMDS = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_FMDS, "src", "models", "Models.jl"))
end

const _FMDS = Models.FreezeoutMesonDensityScan

@testset "FreezeoutMesonDensityScan" begin
    @testset "default output path exists as string" begin
        @test _FMDS.DEFAULT_FREEZEOUT_MESON_DENSITY_OUTPUT_PATH isa String
        @test !isempty(_FMDS.DEFAULT_FREEZEOUT_MESON_DENSITY_OUTPUT_PATH)
    end

    @testset "entrypoint exists" begin
        @test isdefined(_FMDS, :run_freezeout_meson_density_scan)
        @test _FMDS.run_freezeout_meson_density_scan isa Function
        @test Models.run_freezeout_meson_density_scan isa Function
    end

    @testset "header keeps resume-relevant key columns first" begin
        @test startswith(_FMDS.HEADER, "sqrt_s_NN_GeV,muB_MeV,xi,freezeout_profile,path_profile,path_segment,flavor_chemical_profile,meson_chemical_profile,regime,")
    end

    @testset "path-profile columns are present" begin
        @test occursin("path_profile", _FMDS.HEADER)
        @test occursin("path_segment", _FMDS.HEADER)
    end
end
