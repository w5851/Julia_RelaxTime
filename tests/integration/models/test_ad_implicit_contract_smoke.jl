using Test
using StaticArrays
using ForwardDiff

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "models AD implicit contract smoke" begin
    @test isdefined(Models, :build_pnjl_fixedmu_adapters)
    @test isdefined(Models, :build_pnjl_flavor_mu_adapters)

    model = Models.create_model(:PNJL)

    @testset "fixed-mu adapter dual-safe conditions" begin
        adapters = Models.build_pnjl_fixedmu_adapters(model; p_num=12, t_num=4, xi=0.0)
        θ = [0.5, 0.12]
        x, meta = adapters.forward_solve(θ)

        @test length(x) == 5
        @test meta === nothing

        J = ForwardDiff.jacobian(t -> adapters.conditions(t, x, meta), θ)
        @test size(J) == (5, 2)
        @test all(isfinite, J)
    end

    @testset "flavor-mu adapter dual-safe conditions" begin
        adapters = Models.build_pnjl_flavor_mu_adapters(model; p_num=12, t_num=4, xi=0.0)
        θ = [0.5, 0.18, 0.12, 0.06]
        x, meta = adapters.forward_solve(θ)

        @test length(x) == 5
        @test meta === nothing

        J = ForwardDiff.jacobian(t -> adapters.conditions(t, x, meta), θ)
        @test size(J) == (5, 4)
        @test all(isfinite, J)
    end
end
