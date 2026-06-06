using Test

const PROJECT_ROOT_MCP = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_MCP, "src", "models", "Models.jl"))
end

const _MCP = Models.MesonChemicalProfiles

@testset "MesonChemicalProfiles" begin
    @testset "default profile loads" begin
        profile = _MCP.load_meson_chemical_profile()
        @test profile.profile_name == "default"
        @test profile.charge_resolved == false
        @test profile.d_pi == 3
        @test profile.d_K == 4
        @test isapprox(profile.mu_pi_MeV, 0.0; atol=0, rtol=0)
        @test profile.mu_pi_rule === :constant
        @test profile.mu_K_rule === :constant
    end

    @testset "legacy positional constructor defaults pion rule" begin
        profile = _MCP.MesonChemicalProfile(
            "legacy",
            "test",
            "pi",
            "K",
            :pi,
            :K,
            false,
            :constant,
            0.0,
            0.0,
            3,
            4,
        )
        @test profile.mu_pi_rule === :constant
        @test profile.mu_K_rule === :constant
    end

    @testset "charged profile loads" begin
        profile = _MCP.load_meson_chemical_profile(profile="blaschke2019_mu_pi_100")
        @test profile.profile_name == "blaschke2019_mu_pi_100"
        @test profile.charge_resolved == true
        @test profile.pi_label == "pi_minus"
        @test profile.k_label == "K_minus"
        @test profile.d_pi == 1
        @test profile.d_K == 1
        @test isapprox(profile.mu_pi_MeV, 100.0; atol=0, rtol=0)
    end

    @testset "fm conversion helper returns finite values" begin
        profile = _MCP.load_meson_chemical_profile(profile="blaschke2019_mu_pi_134p5")
        out = _MCP.meson_chemical_profile_fm(profile)
        @test out.profile_name == "blaschke2019_mu_pi_134p5"
        @test isfinite(out.mu_pi_fm)
        @test isfinite(out.mu_K_fm)
        @test out.d_pi == 1
    end

    @testset "BU2020 charged signed muK profiles load" begin
        flavor = Models.FlavorChemicalProfiles.flavor_mu_profile_MeV(
            Models.FlavorChemicalProfiles.load_flavor_chemical_profile(profile="bu2020_mu_s_0p2"),
            350.0,
        )
        plus = _MCP.load_meson_chemical_profile(profile="bu2020_kplus_over_piplus_mu_pi_100")
        minus = _MCP.load_meson_chemical_profile(profile="bu2020_kminus_over_piminus_mu_pi_100")
        plus_fm = _MCP.meson_chemical_profile_fm(plus; flavor_mev=flavor)
        minus_fm = _MCP.meson_chemical_profile_fm(minus; flavor_mev=flavor)

        @test plus.pi_channel === :pi_plus
        @test plus.k_channel === :K_plus
        @test minus.pi_channel === :pi_minus
        @test minus.k_channel === :K_minus
        @test plus.mu_K_rule === :mu_u_minus_mu_s_signed
        @test plus_fm.mu_K_fm ≈ 280.0 / Main.Constants_PNJL.ħc_MeV_fm
        @test minus_fm.mu_K_fm ≈ -280.0 / Main.Constants_PNJL.ħc_MeV_fm
    end

    @testset "FixedAsymmetricRho signed pion/kaon profiles resolve flavor mu" begin
        flavor = (mu_u_MeV=270.0, mu_d_MeV=310.0, mu_s_MeV=45.0)
        plus = _MCP.load_meson_chemical_profile(profile="asymmetric_kplus_over_piplus_signed")
        minus = _MCP.load_meson_chemical_profile(profile="asymmetric_kminus_over_piminus_signed")
        plus_fm = _MCP.meson_chemical_profile_fm(plus; flavor_mev=flavor)
        minus_fm = _MCP.meson_chemical_profile_fm(minus; flavor_mev=flavor)

        @test plus.mu_pi_rule === :mu_u_minus_mu_d_signed
        @test plus.mu_K_rule === :mu_u_minus_mu_s_signed
        @test plus_fm.mu_pi_fm ≈ -40.0 / Main.Constants_PNJL.ħc_MeV_fm
        @test minus_fm.mu_pi_fm ≈ 40.0 / Main.Constants_PNJL.ħc_MeV_fm
        @test plus_fm.mu_K_fm ≈ 225.0 / Main.Constants_PNJL.ħc_MeV_fm
        @test minus_fm.mu_K_fm ≈ -225.0 / Main.Constants_PNJL.ħc_MeV_fm
    end

    @testset "signed flavor rules require charged channels" begin
        flavor = (mu_u_MeV=270.0, mu_d_MeV=310.0, mu_s_MeV=45.0)
        neutral_pi = _MCP.MesonChemicalProfile(
            "bad_signed_pi",
            "test",
            "pi",
            "K_plus",
            :pi,
            :K_plus,
            true,
            :mu_u_minus_mu_d_signed,
            :mu_u_minus_mu_s_signed,
            0.0,
            0.0,
            1,
            1,
        )
        neutral_k = _MCP.MesonChemicalProfile(
            "bad_signed_k",
            "test",
            "pi_plus",
            "K",
            :pi_plus,
            :K,
            true,
            :mu_u_minus_mu_d_signed,
            :mu_u_minus_mu_s_signed,
            0.0,
            0.0,
            1,
            1,
        )

        @test_throws ArgumentError _MCP.meson_chemical_profile_fm(neutral_pi; flavor_mev=flavor)
        @test_throws ArgumentError _MCP.meson_chemical_profile_fm(neutral_k; flavor_mev=flavor)
    end
end
