# 体粘滞系数测试
#
# 测试内容：
# 1. bulk_viscosity_coefficients 返回类型检查（无 Dual 泄漏）
# 2. 中间结果验证（与有限差分对比）

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))
using .PNJL: bulk_viscosity_coefficients, solve, FixedMu
using .PNJL.ThermoDerivatives: get_thermal_nodes, set_config, IMPLICIT_SOLVER
using .PNJL.Thermodynamics: calculate_thermo, calculate_rho

using StaticArrays
using ForwardDiff
using LinearAlgebra: dot, norm

@testset "bulk_viscosity_coefficients type check" begin
    T_fm = 150.0 / 197.327
    μ_fm = 300.0 / 197.327

    result = bulk_viscosity_coefficients(T_fm, μ_fm)

    # 检查所有返回值都是 Float64（无 Dual 类型泄漏）
    @test typeof(result.v_n_sq) == Float64
    @test typeof(result.dμB_dT_sigma) == Float64
    @test eltype(result.masses) == Float64
    @test eltype(result.dM_dT) == Float64
    @test eltype(result.dM_dμB) == Float64
    @test typeof(result.s) == Float64
    @test typeof(result.n_B) == Float64

    # 检查数值有限
    @test isfinite(result.v_n_sq)
    @test isfinite(result.dμB_dT_sigma)
    @test all(isfinite.(result.masses))
    @test all(isfinite.(result.dM_dT))
    @test all(isfinite.(result.dM_dμB))
end

@testset "bulk_viscosity_coefficients intermediate verification" begin
    # 使用稳定参数点（远离相变区域）
    # T=150 MeV, μ=300 MeV 在相变附近，有限差分不稳定
    T_fm = 100.0 / 197.327
    μ_fm = 100.0 / 197.327
    ε = 1e-7  # 使用更小的 ε 提高精度

    set_config(xi=0.0, p_num=64, t_num=16)
    thermal_nodes = get_thermal_nodes(64, 16)

    function solve_state(θ_in)
        (x_out, _) = IMPLICIT_SOLVER(θ_in)
        return collect(x_out)
    end

    θ = [T_fm, μ_fm]
    x_base = solve_state(θ)
    x_sv = SVector{5}(Tuple(x_base))
    mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)

    # 计算 dx/dθ
    dx_dθ = ForwardDiff.jacobian(solve_state, θ)

    # 使用有限差分验证 dx/dT
    x_T_plus = solve_state([T_fm + ε, μ_fm])
    x_T_minus = solve_state([T_fm - ε, μ_fm])
    dx_dT_fd = (x_T_plus - x_T_minus) / (2ε)
    @test norm(dx_dθ[:, 1] - dx_dT_fd) / norm(dx_dT_fd) < 0.01

    # 使用有限差分验证 dx/dμ
    x_mu_plus = solve_state([T_fm, μ_fm + ε])
    x_mu_minus = solve_state([T_fm, μ_fm - ε])
    dx_dμ_fd = (x_mu_plus - x_mu_minus) / (2ε)
    @test norm(dx_dθ[:, 2] - dx_dμ_fd) / norm(dx_dμ_fd) < 0.01

    # 计算 ∂s/∂x 和链式法则
    function s_of_x(x_vec)
        x_s = SVector{5}(Tuple(x_vec))
        _, _, s_val, _ = calculate_thermo(x_s, mu_vec, T_fm, thermal_nodes, 0.0)
        return s_val
    end
    ds_dx = ForwardDiff.gradient(s_of_x, x_base)

    function s_of_T(T)
        _, _, s_val, _ = calculate_thermo(x_sv, mu_vec, T, thermal_nodes, 0.0)
        return s_val
    end
    s_T_partial = ForwardDiff.derivative(s_of_T, T_fm)

    ds_dT_total = s_T_partial + dot(ds_dx, dx_dθ[:, 1])

    # 使用有限差分验证总导数
    x_sv_T_plus = SVector{5}(Tuple(x_T_plus))
    x_sv_T_minus = SVector{5}(Tuple(x_T_minus))
    _, _, s_T_plus_total, _ = calculate_thermo(x_sv_T_plus, mu_vec, T_fm + ε, thermal_nodes, 0.0)
    _, _, s_T_minus_total, _ = calculate_thermo(x_sv_T_minus, mu_vec, T_fm - ε, thermal_nodes, 0.0)
    ds_dT_total_fd = (s_T_plus_total - s_T_minus_total) / (2ε)
    
    @test abs(ds_dT_total - ds_dT_total_fd) / abs(ds_dT_total_fd) < 0.01

    # 验证最终结果与函数输出一致
    result = bulk_viscosity_coefficients(T_fm, μ_fm)
    @test isfinite(result.v_n_sq)
    @test isfinite(result.dμB_dT_sigma)
end
