"""
MesonChemicalProfiles 配置解析测试。
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

const _MODELS_PATH_MCP = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
if !isdefined(Main, :Models)
    Base.include(Main, _MODELS_PATH_MCP)
end

using Main.Models.MesonChemicalProfiles: load_meson_chemical_profile, meson_chemical_profile_MeV, meson_chemical_profile_fm
using Main.Models.FlavorChemicalProfiles: load_flavor_chemical_profile, flavor_mu_profile_MeV

@testset "MesonChemicalProfiles charged-channel mapping" begin
    default_profile = load_meson_chemical_profile(profile="default")
    @test default_profile.pi_channel === :pi
    @test default_profile.k_channel === :K

    charged_profile = load_meson_chemical_profile(profile="charge_resolved_neutral")
    @test charged_profile.charge_resolved
    @test charged_profile.pi_channel === :pi_minus
    @test charged_profile.k_channel === :K_minus

    profile_mev = meson_chemical_profile_MeV(charged_profile)
    @test profile_mev.pi_channel === :pi_minus
    @test profile_mev.k_channel === :K_minus
    @test profile_mev.d_pi == 1
    @test profile_mev.d_K == 1

    flavor_profile = load_flavor_chemical_profile(profile="friesen2019_mu_s_0p55")
    flavor_mev = flavor_mu_profile_MeV(flavor_profile, 300.0)

    friesen_plus = load_meson_chemical_profile(profile="friesen2019_kplus_over_piplus_mu_pi_135_rule")
    friesen_minus = load_meson_chemical_profile(profile="friesen2019_kminus_over_piminus_mu_pi_135_rule")

    @test friesen_plus.mu_K_rule === :mu_u_minus_mu_s_signed
    @test friesen_minus.mu_K_rule === :mu_u_minus_mu_s_signed

    friesen_plus_fm = meson_chemical_profile_fm(friesen_plus; flavor_mev=flavor_mev)
    friesen_minus_fm = meson_chemical_profile_fm(friesen_minus; flavor_mev=flavor_mev)
    @test friesen_plus_fm.mu_pi_fm > 0.0
    @test friesen_plus_fm.mu_K_fm ≈ (300.0 - 165.0) / Main.Constants_PNJL.ħc_MeV_fm
    @test friesen_minus_fm.mu_K_fm ≈ -(300.0 - 165.0) / Main.Constants_PNJL.ħc_MeV_fm

    bu_flavor_profile = load_flavor_chemical_profile(profile="bu2020_mu_s_0p2")
    bu_flavor_mev = flavor_mu_profile_MeV(bu_flavor_profile, 350.0)
    bu_plus = load_meson_chemical_profile(profile="bu2020_kplus_over_piplus_mu_pi_100")
    bu_minus = load_meson_chemical_profile(profile="bu2020_kminus_over_piminus_mu_pi_134p5")
    bu_plus_fm = meson_chemical_profile_fm(bu_plus; flavor_mev=bu_flavor_mev)
    bu_minus_fm = meson_chemical_profile_fm(bu_minus; flavor_mev=bu_flavor_mev)

    @test bu_plus.pi_channel === :pi_plus
    @test bu_plus.k_channel === :K_plus
    @test bu_minus.pi_channel === :pi_minus
    @test bu_minus.k_channel === :K_minus
    @test bu_plus_fm.mu_K_fm ≈ 280.0 / Main.Constants_PNJL.ħc_MeV_fm
    @test bu_minus_fm.mu_K_fm ≈ -280.0 / Main.Constants_PNJL.ħc_MeV_fm
end
