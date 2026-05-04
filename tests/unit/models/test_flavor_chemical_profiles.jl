using Test

const PROJECT_ROOT_FCP = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_FCP, "src", "models", "Models.jl"))
end

const _FCP = Models.FlavorChemicalProfiles

@testset "FlavorChemicalProfiles" begin
    @testset "default profile is equal-flavor fixed-mu" begin
        profile = _FCP.load_flavor_chemical_profile()
        vals = _FCP.flavor_mu_profile_MeV(profile, 120.0)
        @test vals.profile_name == "default"
        @test vals.mu_u_MeV == 120.0
        @test vals.mu_d_MeV == 120.0
        @test vals.mu_s_MeV == 120.0
    end

    @testset "blaschke2017 profile lowers strange chemical potential" begin
        profile = _FCP.load_flavor_chemical_profile(profile="blaschke2017_mu_s_0p2")
        vals = _FCP.flavor_mu_profile_MeV(profile, 150.0)
        @test vals.mu_u_MeV == 150.0
        @test vals.mu_d_MeV == 150.0
        @test vals.mu_s_MeV == 30.0
    end
end
