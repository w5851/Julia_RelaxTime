using Test

const PROJECT_ROOT_FREEZEOUT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_FREEZEOUT, "src", "models", "Models.jl"))
end

const _FP = Models.FreezeoutProfiles

@testset "FreezeoutProfiles" begin
    @testset "default profile loads" begin
        profile = _FP.load_freezeout_profile()
        @test profile.profile_name == "default"
        @test profile.family == "cleymans_polynomial_rational_v1"
        @test isapprox(profile.a_GeV, 0.166; atol=0, rtol=0)
        @test isapprox(profile.d_GeV, 1.308; atol=0, rtol=0)
    end

    @testset "named profile loads" begin
        profile = _FP.load_freezeout_profile(profile="cleymans2006")
        @test profile.profile_name == "cleymans2006"
        @test profile.source_tag == "cleymans2006_like"
    end

    @testset "point mapping stays physical" begin
        profile = _FP.load_freezeout_profile()
        pt = _FP.freezeout_point_MeV(profile, 7.7)
        @test pt.sqrt_s_NN_GeV == 7.7
        @test pt.muB_MeV > 0.0
        @test pt.T_MeV > 0.0
        @test pt.T_MeV < 200.0
        @test pt.muq_MeV == pt.muB_MeV / 3.0
    end

    @testset "scan points support traversal for continuation order" begin
        profile = _FP.load_freezeout_profile()
        points = _FP.build_freezeout_scan_points([20.0, 2.5, 7.7]; profile=profile, traversal=:muB_descending)
        @test length(points) == 3
        @test all(diff([p.muB_GeV for p in points]) .< 0.0)

        points_s = _FP.build_freezeout_scan_points([20.0, 2.5, 7.7]; profile=profile, traversal=:sqrts_ascending)
        @test all(diff([p.sqrt_s_NN_GeV for p in points_s]) .> 0.0)
    end

    @testset "invalid traversal rejected" begin
        profile = _FP.load_freezeout_profile()
        @test_throws ArgumentError _FP.build_freezeout_scan_points([7.7, 11.5]; profile=profile, traversal=:bad_order)
    end
end
