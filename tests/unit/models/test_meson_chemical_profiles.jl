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
end
