using Test

# New models entry
_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

using StaticArrays

@testset "PNJLModel.solve_gap models path smoke" begin
    m = Models.create_model(:PNJL; profile="default", physics_profile="default")

    T = 0.9
    mu_vec = @SVector [0.0, 0.0, 0.0]
    xi = 0.0

    p_num = 16
    t_num = 8

    x0 = @SVector [0.02, 0.02, 0.03, 0.5, 0.5]

    solver_models = Models.NLsolveGapSolver(method=:trust_region, jacobian=:finite, xtol=1e-10, ftol=1e-10)

    st_default = Models.solve_gap(m, T, mu_vec;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
    )
    st_custom = Models.solve_gap(m, T, mu_vec;
        solver=solver_models,
        initial_guess=x0,
        residual_norm_max=1e-4,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
    )
    v_default = Models.state_vector(st_default)
    v_custom = Models.state_vector(st_custom)
    @test all(isfinite, v_default)
    @test all(isfinite, v_custom)

    r_default = Models.gap_residual(m, st_default, T, mu_vec; xi=xi, p_num=p_num, t_num=t_num)
    r_custom = Models.gap_residual(m, st_custom, T, mu_vec; xi=xi, p_num=p_num, t_num=t_num)
    @test maximum(abs.(r_default)) <= 1e-2
    @test maximum(abs.(r_custom)) <= 1e-2

    ω_default = Models.omega(m, st_default, T, mu_vec; xi=xi, p_num=p_num, t_num=t_num)
    ω_custom = Models.omega(m, st_custom, T, mu_vec; xi=xi, p_num=p_num, t_num=t_num)

    @test isfinite(ω_default)
    @test isfinite(ω_custom)

    scale = max(1.0, abs(ω_default))
    @test abs(ω_custom - ω_default) / scale <= 1e-2

    # Asymmetric chemical potential is models-only for now; :auto should route to models path.
    mu_asym = @SVector [0.04, 0.00, -0.02]
    st_asym = Models.solve_gap(m, T, mu_asym;
        solver=solver_models,
        initial_guess=x0,
        residual_norm_max=1e-4,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
    )

    v_asym = Models.state_vector(st_asym)
    @test all(isfinite, v_asym)

    r_asym = Models.gap_residual(m, st_asym, T, mu_asym; xi=xi, p_num=p_num, t_num=t_num)
    @test maximum(abs.(r_asym)) <= 1e-2
end
