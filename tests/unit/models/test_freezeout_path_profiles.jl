using Test

const PROJECT_ROOT_FPP = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_FPP, "src", "models", "Models.jl"))
end

const _FPP = Models.FreezeoutPathProfiles

@testset "FreezeoutPathProfiles" begin
    @testset "baseline profile preserves freezeout path" begin
        freezeout_profile = Models.load_freezeout_profile()
        path_profile = _FPP.load_freezeout_path_profile()
        points = _FPP.build_freezeout_path_points([7.7, 20.0];
            freezeout_profile=freezeout_profile,
            path_profile=path_profile,
            traversal=:sqrts_ascending,
        )
        @test length(points) == 2
        @test all(p -> p.path_profile == "baseline_freezeout", points)
        @test all(p -> p.path_segment == "baseline_freezeout", points)
    end

    @testset "const-T proxy switches high-energy segment" begin
        freezeout_profile = Models.load_freezeout_profile()
        path_profile = _FPP.load_freezeout_path_profile(profile="freezeout_plus_constT_proxy")
        points = _FPP.build_freezeout_path_points([7.7, 20.0];
            freezeout_profile=freezeout_profile,
            path_profile=path_profile,
            traversal=:sqrts_ascending,
        )
        @test points[1].path_segment == "baseline_freezeout"
        @test points[2].path_segment == "constant_T"
        @test points[2].T_MeV == 166.0
    end

    @testset "critical-line-plus-constT proxy is available as named profile" begin
        profile = _FPP.load_freezeout_path_profile(profile="critical_line_plus_constT_proxy")
        @test profile.profile_name == "critical_line_plus_constT_proxy"
        @test profile.family == "freezeout_plus_constT_v1"
        @test profile.switch_sqrt_s_NN_GeV == 8.0
        @test profile.constant_T_MeV == 166.0
    end
end
