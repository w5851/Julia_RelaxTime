using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Models implicit differentiation flavor-mu smoke" begin
    points = (
        (T_fm=0.48, mu_vec=SVector(0.0, 0.0, 0.0)),
        (T_fm=0.52, mu_vec=SVector(0.18, 0.12, 0.06)),
        (T_fm=0.64, mu_vec=SVector(0.24, 0.10, 0.02)),
    )

    for point in points
        result = Models.solve_pnjl_with_flavor_mu_derivatives(
            point.T_fm,
            point.mu_vec;
            p_num=12,
            t_num=4,
            xi=0.0,
        )

        @test length(result.x) == 5
        @test result.mu_vec == point.mu_vec
        @test length(result.dx_dT) == 5
        @test size(result.dx_dmu_vec) == (5, 3)
        @test all(isfinite, result.x)
        @test all(isfinite, result.dx_dT)
        @test all(isfinite, result.dx_dmu_vec)
        @test all(isfinite, Models.symmetric_mu_direction_derivative(result.dx_dmu_vec))
    end

    symmetric_mu = 0.14
    symmetric_result = Models.solve_pnjl_with_flavor_mu_derivatives(
        0.5,
        SVector(symmetric_mu, symmetric_mu, symmetric_mu);
        p_num=12,
        t_num=4,
        xi=0.0,
    )
    scalar_result = Models.solve_pnjl_with_derivatives(
        0.5,
        symmetric_mu;
        p_num=12,
        t_num=4,
        xi=0.0,
    )

    @test all(isapprox.(symmetric_result.x, scalar_result.x; rtol=1e-6, atol=1e-8))
    @test all(isapprox.(symmetric_result.dx_dT, scalar_result.dx_dT; rtol=1e-6, atol=1e-8))
    @test all(isapprox.(Models.symmetric_mu_direction_derivative(symmetric_result.dx_dmu_vec), scalar_result.dx_dμ; rtol=1e-6, atol=1e-8))
end