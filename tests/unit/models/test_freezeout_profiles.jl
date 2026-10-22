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

    @testset "published default chemical freezeout lines" begin
        cleymans = _FP.load_freezeout_profile(profile="cleymans2006")
        andronic = _FP.load_freezeout_profile(profile="andronic2018")
        @test andronic isa _FP.AndronicFreezeoutProfile
        @test cleymans isa _FP.FreezeoutParameterProfile
        @test andronic isa _FP.AbstractFreezeoutParameterProfile
        # Independently evaluated published formulas, in MeV (rounded here).
        for (s, Tc, muc, Ta, mua) in (
            (7.7, 139.612, 421.650, 138.428, 406.359),
            (11.5, 151.593, 315.980, 149.552, 303.224),
            (19.6, 160.008, 205.958, 155.585, 196.770),
            (200.0, 165.923, 23.525, 158.384, 22.312),
        )
            pc = _FP.freezeout_point_MeV(cleymans, s)
            pa = _FP.freezeout_point_MeV(andronic, s)
            @test pc.T_MeV ≈ Tc atol=0.001
            @test pc.muB_MeV ≈ muc atol=0.001
            @test pa.T_MeV ≈ Ta atol=0.001
            @test pa.muB_MeV ≈ mua atol=0.001
            @test _FP.freezeout_temperature_GeV(andronic, pa.muB_GeV) ≈ pa.T_GeV rtol=1e-14
            @test pa.T_fm * Main.Constants_PNJL.ħc_MeV_fm ≈ pa.T_MeV
            @test pa.muB_fm * Main.Constants_PNJL.ħc_MeV_fm ≈ pa.muB_MeV
            # Legacy default mapping must remain bit-for-bit identical.
            old_mu = 1.308 / (1 + 0.273s)
            old_T = 0.166 - 0.139old_mu^2 - 0.053old_mu^4
            @test _FP.freezeout_point_GeV(_FP.load_freezeout_profile(), s) ==
                (sqrt_s_NN_GeV=s, muB_GeV=old_mu, T_GeV=old_T)
        end
        @test _FP.freezeout_temperature_GeV(andronic, 0.) == 0.1584
        @test _FP.freezeout_temperature_GeV(andronic, andronic.d_GeV) == 0.
        @test _FP.freezeout_point_GeV(andronic, 0.).T_GeV == 0.
        @test _FP.freezeout_point_GeV(andronic, 1e200).T_GeV == 0.1584
        path = Models.FreezeoutPathProfiles.build_freezeout_path_points([200., 7.7];
            freezeout_profile=andronic, traversal=:sqrts_ascending)
        @test first(path).T_MeV ≈ 138.428 atol=0.001
        @test first(path).path_segment == "baseline_freezeout"
        for profile in (cleymans, andronic), invalid in (-1., NaN, Inf)
            @test_throws ArgumentError _FP.freezeout_point_MeV(profile, invalid)
            @test_throws ArgumentError _FP.freezeout_temperature_GeV(profile, invalid)
        end
        @test_throws ArgumentError _FP.freezeout_temperature_GeV(andronic, 1.4)
    end

    @testset "bad named profiles cannot silently select another line" begin
        @test_throws ArgumentError _FP.load_freezeout_profile(profile="not_a_freezeout_profile")
        mktempdir() do dir
            write(joinpath(dir, "bad.toml"), "invalid = [")
            @test_throws ArgumentError _FP.load_freezeout_profile(profile="bad", dir=dir)
            original = read(joinpath(_FP.freezeout_profile_dir(), "andronic2018.toml"), String)
            for changed in (
                replace(original, "andronic_logistic_rational_v1" => "unknown_family"),
                replace(original, "d_GeV = 1.3075" => ""),
                replace(original, "log_width = 0.45" => "log_width = 0.0"),
                replace(original, "T_lim_GeV = 0.1584" => "T_lim_GeV = nan"),
            )
                write(joinpath(dir, "bad.toml"), changed)
                @test_throws ArgumentError _FP.load_freezeout_profile(profile="bad", dir=dir)
            end
        end
    end
end
