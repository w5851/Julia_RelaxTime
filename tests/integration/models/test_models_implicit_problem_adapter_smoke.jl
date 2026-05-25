# Models implicit problem adapter smoke test
#
# Validates that Models-side implicit residual adapters remain usable without
# constructing retired implicit differentiation factories.

using Test
using ForwardDiff
using LinearAlgebra

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

_models_entry = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    Base.include(Main, _models_entry)
end

@testset "Models implicit problem adapter (NJL)" begin
    model = Models.create_model(:NJL)

    problem = Models.build_njl_problem(model; p_num=10, t_num=4, xi=0.0)

    θ = [0.12, 0.0]
    x, meta = problem.forward_solve(θ)
    residual = problem.conditions(θ, x, meta)

    @test length(x) == 3
    @test all(isfinite, x)
    @test length(residual) == 3
    @test all(isfinite, residual)
    @test norm(residual) <= 1e-5

    J = ForwardDiff.jacobian(θv -> problem.conditions(θv, x, meta), θ)
    @test size(J) == (3, 2)
    @test all(isfinite, J)
end
