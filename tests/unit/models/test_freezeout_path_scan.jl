using Test

const PROJECT_ROOT_FREEZEOUT_SCAN = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_FREEZEOUT_SCAN, "src", "models", "Models.jl"))
end

const _FPSCAN = Models.FreezeoutPathScan

@testset "FreezeoutPathScan" begin
    @testset "default output path exists as string" begin
        @test _FPSCAN.DEFAULT_FREEZEOUT_OUTPUT_PATH isa String
        @test !isempty(_FPSCAN.DEFAULT_FREEZEOUT_OUTPUT_PATH)
    end

    @testset "entrypoint exists" begin
        @test isdefined(_FPSCAN, :run_freezeout_fixedmu_scan)
        @test _FPSCAN.run_freezeout_fixedmu_scan isa Function
        @test Models.run_freezeout_fixedmu_scan isa Function
    end

    @testset "schedule points align with continuation traversal" begin
        profile = Models.load_freezeout_profile()
        points = Models.build_freezeout_scan_points([20.0, 7.7, 2.5]; profile=profile, traversal=:muB_descending)
        @test length(points) == 3
        @test all(diff([p.muB_MeV for p in points]) .< 0.0)
    end

    @testset "header key order keeps resume-compatible columns" begin
        @test startswith(_FPSCAN.HEADER, "sqrt_s_NN_GeV,muB_MeV,xi,")
    end

    @testset "path-profile columns are present" begin
        @test occursin("path_profile", _FPSCAN.HEADER)
        @test occursin("path_segment", _FPSCAN.HEADER)
    end
end
