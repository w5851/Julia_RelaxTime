using Test
using ForwardDiff

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

const _MODELS_ENTRY = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
if !(isdefined(Main, :Models) && isdefined(Main.Models, :create_implicit_gap_solver))
    Base.include(Main, _MODELS_ENTRY)
end

@testset "Models implicit differentiation (NJL2)" begin
    model = Models.create_model(:NJL2)

    f = Models.create_implicit_gap_solver(model; p_num=10, t_num=4, xi=0.0)

    θ = [0.12, 0.0]
    x, _ = f(θ)
    @test length(x) == 2
    @test all(isfinite, x)

    J = ForwardDiff.jacobian(θv -> f(θv)[1], θ)
    @test size(J) == (2, 2)
    @test all(isfinite, J)
end
