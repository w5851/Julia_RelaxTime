# PNJL TaylorDiff derivative Jacobian smoke test
#
# The old ImplicitDifferentiation.jl integration path has been retired. This
# file now checks the same implicit-theorem identity against the TD derivative
# wrapper without constructing an ImplicitFunction.

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()
const PNJL = Models.pnjl_module()
const solve_with_derivatives = getproperty(PNJL, :solve_with_derivatives)
const gap_conditions = getproperty(PNJL, :gap_conditions)
const GapParams = getproperty(PNJL, :GapParams)
const cached_nodes = getproperty(PNJL, :cached_nodes)

using ForwardDiff
using StaticArrays
using LinearAlgebra

@testset "manual implicit theorem basic behavior" begin
    # x[1]^3 - theta[1] = 0 and x[2]^2 - theta[2] = 0 at theta=(8, 4).
    θ = [8.0, 4.0]
    x = [cbrt(θ[1]), sqrt(θ[2])]

    conditions(theta, state) = [state[1]^3 - theta[1], state[2]^2 - theta[2]]
    dF_dx = ForwardDiff.jacobian(state -> conditions(θ, state), x)
    dF_dθ = ForwardDiff.jacobian(theta -> conditions(theta, x), θ)
    dx_dθ = -dF_dx \ dF_dθ

    @test size(dx_dθ) == (2, 2)
    @test isapprox(dx_dθ[1, 1], 1 / (3x[1]^2); rtol=1e-12, atol=1e-12)
    @test isapprox(dx_dθ[2, 2], 1 / (2x[2]); rtol=1e-12, atol=1e-12)
    @test isapprox(norm(dF_dx * dx_dθ + dF_dθ), 0.0; atol=1e-12)
end

@testset "PNJL TD derivative wrapper satisfies implicit residual identity" begin
    T_fm = 0.55
    μ_fm = 0.30
    p_num = 8
    t_num = 4
    xi = 0.0

    d = solve_with_derivatives(
        T_fm,
        μ_fm;
        order=1,
        p_num=p_num,
        t_num=t_num,
        xi=xi,
        series_residual_tol=1e-7,
    )

    thermal_nodes = cached_nodes(p_num, t_num)
    mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)
    x_star = SVector{5}(Tuple(d.x))
    params = GapParams(T_fm, thermal_nodes, xi; p_num=p_num, t_num=t_num, model_kind=:PNJL)
    F_star = gap_conditions(x_star, mu_vec, params)

    @test norm(F_star) <= 1e-6

    F_of_x = x_vec -> Vector(gap_conditions(SVector{5}(Tuple(x_vec)), mu_vec, params))
    dF_dx = ForwardDiff.jacobian(F_of_x, collect(x_star))

    F_of_T = T -> begin
        p = GapParams(T, thermal_nodes, xi; p_num=p_num, t_num=t_num, model_kind=:PNJL)
        Vector(gap_conditions(x_star, mu_vec, p))
    end
    dF_dT = ForwardDiff.derivative(F_of_T, T_fm)

    F_of_mu = μ -> begin
        μv = SVector{3}(μ, μ, μ)
        Vector(gap_conditions(x_star, μv, params))
    end
    dF_dμ = ForwardDiff.derivative(F_of_mu, μ_fm)

    residual_T = dF_dx * d.dx_dT + dF_dT
    residual_μ = dF_dx * d.dx_dμ + dF_dμ

    @test norm(residual_T) <= 1e-4 * (1 + norm(dF_dT))
    @test norm(residual_μ) <= 1e-4 * (1 + norm(dF_dμ))
end
