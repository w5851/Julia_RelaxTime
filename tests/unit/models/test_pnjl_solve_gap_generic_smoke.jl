using Test

# New models entry
_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

using StaticArrays

@testset "Models PNJL generic solve_gap smoke (5D)" begin
    m = Models.create_model(:PNJL; profile="default", physics_profile="default")

    T = 0.9
    mu_vec = @SVector [0.0, 0.0, 0.0]
    xi = 0.0

    p_num = 16
    t_num = 8

    # A mild seed near the omega-bridge smoke state.
    x0 = @SVector [0.02, 0.02, 0.03, 0.5, 0.5]

    # Force calling the generic AbstractPNJLModel method even if PNJLModel has a more specific solver.
    solver = Models.NLsolveGapSolver(method=:trust_region, jacobian=:finite, xtol=1e-10, ftol=1e-10)

    st = Models.solve_gap(
        m,
        T,
        mu_vec;
        solver_backend=:models,
        solver=solver,
        initial_guess=x0,
        residual_norm_max=1e-4,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
    )

    v = Models.state_vector(st)
    @test length(v) == 5
    @test all(isfinite, v)

    r = Models.gap_residual(m, st, T, mu_vec; xi=xi, p_num=p_num, t_num=t_num)
    @test length(r) == 5
    @test all(isfinite, r)

    # Loose acceptance: smoke-level check that it's near a stationary point.
    @test maximum(abs.(r)) <= 1e-2
end
