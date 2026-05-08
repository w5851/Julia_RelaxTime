using Test

const PROJECT_ROOT_IPP = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_IPP, "src", "models", "Models.jl"))
end

const _IPP = Models.IsentropicPathProfiles

@testset "IsentropicPathProfiles" begin
    @testset "default profile binds sigma target" begin
        profile = _IPP.load_isentropic_path_profile(sigma_target=30.0)
        @test profile.profile_name == "default"
        @test profile.family == "fixed_sigma_v1"
        @test profile.sigma_target == 30.0
    end

    @testset "build points keeps sigma metadata and sorting" begin
        points = _IPP.build_isentropic_path_points([170.0, 150.0];
            sigma_target=25.0,
            traversal=:T_ascending,
        )
        @test length(points) == 2
        @test points[1].T_MeV == 150.0
        @test all(pt -> pt.path_family == "isentropic", points)
        @test all(pt -> pt.path_segment == "fixed_sigma", points)
        @test all(pt -> pt.sigma_target == 25.0, points)
        @test points[1].path_point_index == 0
        @test points[2].path_point_index == 1
    end
end
