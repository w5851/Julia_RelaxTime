using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "solver config isolation smoke" begin
    @testset "pnjl implicit solver instances keep independent config" begin
        θ = [0.52, 0.16]

        solver_a = Models.create_pnjl_implicit_solver(xi=0.0, p_num=12, t_num=4)
        solver_b = Models.create_pnjl_implicit_solver(xi=0.35, p_num=24, t_num=8)

        x_a, _ = solver_a(θ)
        x_b, _ = solver_b(θ)

        @test length(x_a) == 5
        @test length(x_b) == 5
        @test all(isfinite, x_a)
        @test all(isfinite, x_b)

        @test !all(isapprox.(x_a, x_b; rtol=1e-10, atol=1e-12))
    end
end
