# Models implicit problem adapter smoke test (NJL2)

using Test
using ForwardDiff
using LinearAlgebra

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

const _MODELS_ENTRY = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    Base.include(Main, _MODELS_ENTRY)
end

@testset "Models implicit problem adapter (NJL2)" begin
    model = Models.create_model(:NJL2)

    problem = Models.build_njl_problem(model; p_num=10, t_num=4, xi=0.0)

    θ = [0.12, 0.0]
    x, meta = problem.forward_solve(θ)
    residual = problem.conditions(θ, x, meta)

    @test length(x) == 2
    @test all(isfinite, x)
    @test length(residual) == 2
    @test all(isfinite, residual)
    @test norm(residual) <= 1e-5

    J = ForwardDiff.jacobian(θv -> problem.conditions(θv, x, meta), θ)
    @test size(J) == (2, 2)
    @test all(isfinite, J)
end
