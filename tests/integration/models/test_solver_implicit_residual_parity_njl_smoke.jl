using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "solver/implicit residual parity njl smoke" begin
    p_num = 10
    t_num = 4
    xi = 0.0

    @testset "NJL2 parity" begin
        model = Models.create_model(:NJL2)
        problem = Models.build_njl_problem(model; p_num=p_num, t_num=t_num, xi=xi)

        theta = [0.12, 0.0]
        x, meta = problem.forward_solve(theta)
        @test meta === nothing
        @test length(x) == 2

        params = Models.GapParams(theta[1], Models.cached_nodes(p_num, t_num), xi;
            p_num=p_num, t_num=t_num, model_kind=:NJL2)
        mu_vec = SVector{3}(theta[2], theta[2], theta[2])

        F_solver = zeros(2)
        Models.gap_core_residual!(F_solver, x, mu_vec, params)

        F_implicit = problem.conditions(theta, x, meta)
        @test length(F_implicit) == 2
        @test all(isfinite, F_implicit)

        for i in 1:2
            @test isapprox(F_solver[i], F_implicit[i]; rtol=1e-8, atol=1e-9)
        end
    end

    @testset "NJL parity" begin
        model = Models.create_model(:NJL)
        problem = Models.build_njl_problem(model; p_num=p_num, t_num=t_num, xi=xi)

        theta = [0.12, 0.0]
        x, meta = problem.forward_solve(theta)
        @test meta === nothing
        @test length(x) == 3

        params = Models.GapParams(theta[1], Models.cached_nodes(p_num, t_num), xi;
            p_num=p_num, t_num=t_num, model_kind=:NJL)
        mu_vec = SVector{3}(theta[2], theta[2], theta[2])

        F_solver = zeros(3)
        Models.gap_core_residual!(F_solver, x, mu_vec, params)

        F_implicit = problem.conditions(theta, x, meta)
        @test length(F_implicit) == 3
        @test all(isfinite, F_implicit)

        for i in 1:3
            @test isapprox(F_solver[i], F_implicit[i]; rtol=1e-8, atol=1e-9)
        end
    end
end
