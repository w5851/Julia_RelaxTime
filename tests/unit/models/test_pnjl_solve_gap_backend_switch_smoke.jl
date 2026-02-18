using Test

# New models entry
_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

using StaticArrays

@testset "PNJLModel.solve_gap solver_backend switch smoke" begin
    m = Models.create_model(:PNJL; profile="default", physics_profile="default")

    T = 0.9
    mu_vec = @SVector [0.0, 0.0, 0.0]
    xi = 0.0

    p_num = 16
    t_num = 8

    x0 = @SVector [0.02, 0.02, 0.03, 0.5, 0.5]

    solver_models = Models.NLsolveGapSolver(method=:trust_region, jacobian=:finite, xtol=1e-10, ftol=1e-10)

    st_legacy = Models.solve_gap(m, T, mu_vec; solver_backend=:legacy, xi=xi, p_num=p_num, t_num=t_num)
    st_models = Models.solve_gap(m, T, mu_vec;
        solver_backend=:models,
        solver=solver_models,
        initial_guess=x0,
        residual_norm_max=1e-4,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
    )

    vL = Models.state_vector(st_legacy)
    vM = Models.state_vector(st_models)
    @test all(isfinite, vL)
    @test all(isfinite, vM)

    rL = Models.gap_residual(m, st_legacy, T, mu_vec; xi=xi, p_num=p_num, t_num=t_num)
    rM = Models.gap_residual(m, st_models, T, mu_vec; xi=xi, p_num=p_num, t_num=t_num)
    @test maximum(abs.(rL)) <= 1e-2
    @test maximum(abs.(rM)) <= 1e-2

    ωL = Models.omega(m, st_legacy, T, mu_vec; xi=xi, p_num=p_num, t_num=t_num)
    ωM = Models.omega(m, st_models, T, mu_vec; xi=xi, p_num=p_num, t_num=t_num)

    @test isfinite(ωL)
    @test isfinite(ωM)

    # Loose comparability check (smoke): energies should be close.
    scale = max(1.0, abs(ωL))
    @test abs(ωM - ωL) / scale <= 1e-2
end
