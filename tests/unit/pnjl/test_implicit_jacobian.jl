# ImplicitDifferentiation.jl 集成测试
#
# 测试内容：
# 1. ImplicitFunction + ForwardDiff.jacobian 的基本行为
# 2. 验证 dx/dθ 是否符合隐函数定理

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
Models.legacy_pnjl_module()
const PNJL = Models.legacy_pnjl_module()
const solve = getproperty(PNJL, :solve)
const FixedMu = getproperty(PNJL, :FixedMu)
const gap_conditions = getproperty(PNJL, :gap_conditions)
const GapParams = getproperty(PNJL, :GapParams)
const cached_nodes = getproperty(PNJL, :cached_nodes)
const ThermoDerivatives = getproperty(PNJL, :ThermoDerivatives)
const IMPLICIT_SOLVER = getproperty(ThermoDerivatives, :IMPLICIT_SOLVER)
const set_config = getproperty(ThermoDerivatives, :set_config)

using ForwardDiff
using NLsolve
using ImplicitDifferentiation
using StaticArrays
using LinearAlgebra

@testset "ImplicitFunction basic behavior" begin
    # 简单多维隐函数问题（解析 forward solve，避免数值求解引入的不稳定/奇异）：
    # x[1]^3 - θ[1] = 0  => x[1] = cbrt(θ[1])
    # x[2]^2 - θ[2] = 0  => x[2] = sqrt(θ[2])
    function forward_solve_multi(θ::AbstractVector)
        x1 = cbrt(θ[1])
        x2 = sqrt(θ[2])
        return ([x1, x2], nothing)
    end

    function conditions_multi(θ::AbstractVector, x::AbstractVector, _)
        return [x[1]^3 - θ[1], x[2]^2 - θ[2]]
    end

    implicit_multi = ImplicitFunction(
        forward_solve_multi, conditions_multi;
        linear_solver=DirectLinearSolver(),
        representation=MatrixRepresentation()
    )

    θ = [8.0, 4.0]
    (x, _) = implicit_multi(θ)
    
    @test length(x) == 2
    @test x[1]^3 ≈ θ[1] atol=1e-8
    @test x[2]^2 ≈ θ[2] atol=1e-8

    # 测试 Jacobian 计算
    function my_func_multi(θ_in)
        (x_out, _) = implicit_multi(θ_in)
        return [x_out[1] + x_out[2], x_out[1]^2, θ_in[1] + θ_in[2]]
    end

    J = ForwardDiff.jacobian(my_func_multi, θ)
    @test size(J) == (3, 2)
    @test all(isfinite.(J))
    @test eltype(J) == Float64
end

@testset "PNJL implicit solver dx/dθ verification" begin
    T_fm = 150.0 / 197.327
    μ_fm = 300.0 / 197.327

    # 求解
    thermal_nodes = cached_nodes(64, 16)
    result = solve(FixedMu(), T_fm, μ_fm)
    x_star = result.x_state
    
    @test result.converged

    # 验证 F(x*, θ) ≈ 0
    mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)
    x_sv = SVector{5}(Tuple(x_star))
    params = GapParams(T_fm, thermal_nodes, 0.0)
    F_star = gap_conditions(x_sv, mu_vec, params)
    
    @test norm(F_star) < 1e-8

    # 计算 ∂F/∂x (Jacobian)
    function F_of_x(x_vec)
        x_s = SVector{5}(Tuple(x_vec))
        return Vector(gap_conditions(x_s, mu_vec, params))
    end
    dF_dx = ForwardDiff.jacobian(F_of_x, collect(x_star))
    
    @test size(dF_dx) == (5, 5)
    @test all(isfinite.(dF_dx))

    # 计算 ∂F/∂T 和 ∂F/∂μ
    function F_of_T(T)
        p = GapParams(T, thermal_nodes, 0.0)
        return Vector(gap_conditions(x_sv, mu_vec, p))
    end
    dF_dT = ForwardDiff.derivative(F_of_T, T_fm)

    function F_of_mu(μ)
        mu_v = SVector{3}(μ, μ, μ)
        return Vector(gap_conditions(x_sv, mu_v, params))
    end
    dF_dμ = ForwardDiff.derivative(F_of_mu, μ_fm)

    # 隐函数定理：dx/dθ = -(∂F/∂x)⁻¹ · (∂F/∂θ)
    dx_dT_implicit = -dF_dx \ dF_dT
    dx_dμ_implicit = -dF_dx \ dF_dμ

    # 与 ImplicitDifferentiation.jl 的结果比较
    set_config(xi=0.0, p_num=64, t_num=16)

    function solve_state(θ_in)
        (x_out, _) = IMPLICIT_SOLVER(θ_in)
        return collect(x_out)
    end

    θ = [T_fm, μ_fm]
    dx_dθ_lib = ForwardDiff.jacobian(solve_state, θ)

    # 验证相对误差 < 1%
    @test norm(dx_dT_implicit - dx_dθ_lib[:, 1]) / norm(dx_dT_implicit) < 0.01
    @test norm(dx_dμ_implicit - dx_dθ_lib[:, 2]) / norm(dx_dμ_implicit) < 0.01
end
