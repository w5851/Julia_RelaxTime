# Models implicit differentiation smoke test
#
# Validates that Models-side ImplicitDifferentiation wiring works:
# - create_implicit_gap_solver(::NJLModel)
# - primal solve (θ -> x)
# - dx/dθ via ForwardDiff.jacobian (implicit differentiation)

using Test
using ForwardDiff

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

_models_entry = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
if !(isdefined(Main, :Models) && isdefined(Main.Models, :create_implicit_gap_solver))
    Base.include(Main, _models_entry)
end

@testset "Models implicit differentiation (NJL)" begin
    model = Models.create_model(:NJL)

    # Keep this intentionally cheap; just validates the plumbing.
    f = Models.create_implicit_gap_solver(model; p_num=10, t_num=4, xi=0.0)

    θ = [0.12, 0.0]
    x, _ = f(θ)
    @test length(x) == 3
    @test all(isfinite, x)

    J = ForwardDiff.jacobian(θv -> f(θv)[1], θ)
    @test size(J) == (3, 2)
    @test all(isfinite, J)
end
