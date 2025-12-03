# 保证 GaussLegendre 模块已加载，供积分节点生成复用
const _GAUSSLEGENDRE_PATH = normpath(joinpath(@__DIR__, "..", "..", "integration", "GaussLegendre.jl"))
if !isdefined(Main, :GaussLegendre)
    Base.include(Main, _GAUSSLEGENDRE_PATH)
end

module AnisoGapSolver

using LinearAlgebra
using Base.MathConstants: π
using StaticArrays
using ForwardDiff
using NLsolve
using SpecialFunctions: log
using Main.GaussLegendre:
    gauleg,
    DEFAULT_COSΘ_HALF_NODES,
    DEFAULT_COSΘ_HALF_WEIGHTS,
    DEFAULT_MOMENTUM_NODES,
    DEFAULT_MOMENTUM_WEIGHTS

using ..Constants_PNJL:
    ħc_MeV_fm,
    α,
    N_color,
    N_flavor,
    ρ0_inv_fm3,
    m_ud0_inv_fm,
    m_s0_inv_fm,
    Λ_inv_fm,
    G_fm2,
    K_fm5,
    T0_inv_fm,
    a0,
    a1,
    a2,
    b3,
    b4

const DEFAULT_THETA_COUNT = length(DEFAULT_COSΘ_HALF_NODES)
const DEFAULT_MOMENTUM_COUNT = length(DEFAULT_MOMENTUM_NODES)
const THETA_DEFAULT_NODES = DEFAULT_COSΘ_HALF_NODES
const THETA_DEFAULT_WEIGHTS = DEFAULT_COSΘ_HALF_WEIGHTS .* 2.0
const THERMAL_DEFAULT_NODES = DEFAULT_MOMENTUM_NODES
const THERMAL_DEFAULT_WEIGHTS = DEFAULT_MOMENTUM_WEIGHTS

const ρ0 = ρ0_inv_fm3
const DEFAULT_RHO_GUESS = [-0.07441, -0.07441, -1.86717, 0.02032, 0.02372, 2.15664, 2.15664, 2.15664]
const DEFAULT_MU_GUESS = [-1.8, -1.8, -2.1, 0.8, 0.8]
const NODE_CACHE = Dict{Tuple{Int, Int}, NTuple{3, Matrix{Float64}}}()

export SolverResult, solve_fixed_rho, solve_fixed_mu, DEFAULT_RHO_GUESS, DEFAULT_MU_GUESS

struct SolverResult
    mode::Symbol
    converged::Bool
    solution::Vector{Float64}
    mu_vec::SVector{3, Float64}
    pressure::Float64
    rho::Float64
    entropy::Float64
    energy::Float64
    iterations::Int
    residual_norm::Float64
    xi::Float64
end

@inline function safe_log(x; min_val=1e-16)
    x <= 0 && return log(min_val)
    return x < min_val ? log(min_val) : log(x)
end

@inline function theta_nodes(t_num::Int)
    if t_num == DEFAULT_THETA_COUNT
        return THETA_DEFAULT_NODES, THETA_DEFAULT_WEIGHTS
    end
    nodes, weights = gauleg(0.0, 1.0, t_num)
    return nodes, weights .* 2.0
end

@inline function thermal_nodes(p_num::Int)
    if p_num == DEFAULT_MOMENTUM_COUNT
        return THERMAL_DEFAULT_NODES, THERMAL_DEFAULT_WEIGHTS
    end
    return gauleg(0.0, 10.0, p_num)
end

function build_nodes(p_num::Int, t_num::Int)
    momentum_nodes, momentum_weights = thermal_nodes(p_num)
    cosθ_nodes, cosθ_weights = theta_nodes(t_num)

    thermal_p_mesh = repeat(momentum_nodes, 1, t_num)
    cosθ_mesh = repeat(cosθ_nodes', p_num, 1)
    weight_mesh = momentum_weights * cosθ_weights'
    thermal_coefficients = weight_mesh .* thermal_p_mesh.^2 ./ (2 * π)^2

    return (thermal_p_mesh, cosθ_mesh, thermal_coefficients)
end

function cached_nodes(p_num::Int, t_num::Int)
    key = (p_num, t_num)
    return get!(NODE_CACHE, key) do
        build_nodes(p_num, t_num)
    end
end

@inline function calculate_chiral(φ::SVector{3, TF}) where {TF}
    2 * G_fm2 * sum(φ .^ 2) - 4 * K_fm5 * prod(φ)
end

@inline function calculate_U(T_fm::Real, Φ::Real, Φ̄::Real)
    T_ratio = T0_inv_fm / T_fm
    Ta = a0 + a1 * T_ratio + a2 * T_ratio^2
    Tb = b3 * T_ratio^3
    value = 1 - 6 * Φ̄ * Φ + 4 * (Φ̄^3 + Φ^3) - 3 * (Φ̄ * Φ)^2
    return T_fm^4 * (-0.5 * Ta * Φ̄ * Φ + Tb * safe_log(value))
end

@inline function calculate_mass_vec(φ::SVector{3, TF}) where {TF}
    φ_u, φ_d, φ_s = φ
    return SVector{3, TF}(
        m_ud0_inv_fm - 4 * G_fm2 * φ_u + 2 * K_fm5 * φ_d * φ_s,
        m_ud0_inv_fm - 4 * G_fm2 * φ_d + 2 * K_fm5 * φ_u * φ_s,
        m_s0_inv_fm - 4 * G_fm2 * φ_s + 2 * K_fm5 * φ_u * φ_d,
    )
end

@inline function calculate_energy_isotropic(mass_i, p)
    return sqrt(p^2 + mass_i^2)
end

@inline function calculate_energy_anisotropic(mass_i, p, xi, t)
    return sqrt(p^2 + mass_i^2 + xi * (p * t)^2)
end

@inline function calculate_log_term(E_i, mu_i, T_fm, Φ, Φ̄)
    invT = 1.0 / T_fm
    x_i = (E_i - mu_i) * invT
    x_i_anti = (E_i + mu_i) * invT

    exp1 = exp(-x_i)
    exp2 = exp1 * exp1
    exp3 = exp1 * exp2
    exp1_anti = exp(-x_i_anti)
    exp2_anti = exp1_anti * exp1_anti
    exp3_anti = exp1_anti * exp2_anti

    f1_val = 1.0 + 3.0 * Φ * exp1 + 3.0 * Φ̄ * exp2 + exp3
    f2_val = 1.0 + 3.0 * Φ̄ * exp1_anti + 3.0 * Φ * exp2_anti + exp3_anti

    return safe_log(f1_val) + safe_log(f2_val)
end

"""计算真空项的解析积分, 对应 -2Nc ∑∫ d^3p/(2π)^3 E_i"""
@inline function vacuum_integral(mass::TF) where {TF}
    Λ = convert(TF, Λ_inv_fm)
    mass_abs = abs(mass)
    epsilon = one(mass_abs) * 1e-12
    mass_safe = mass_abs + epsilon
    sqrt_term = sqrt(Λ^2 + mass_safe^2)
    poly_part = Λ * sqrt_term * (2 * Λ^2 + mass_safe^2)
    log_term = mass_safe^4 * log((Λ + sqrt_term) / mass_safe)
    return (poly_part - log_term) / (16 * π^2)
end

function calculate_energy_sum(masses::SVector{3, TF}) where {TF}
    total = zero(TF)
    @inbounds for i in 1:3
        total += vacuum_integral(masses[i])
    end
    return -2 * N_color * total
end

"""计算RS动量各向异性分布下巨热力学势热项(对数和)"""
function calculate_log_sum(masses::SVector{3, TF}, p_nodes, cosθ_nodes, coefficients, Φ, Φ̄, mu_vec, T_fm, xi) where {TF}
    total = zero(TF)
    @inbounds for i in 1:3
        mass_i = masses[i]
        mu_i = mu_vec[i]
        for idx in eachindex(p_nodes)
            E_i = calculate_energy_anisotropic(mass_i, p_nodes[idx], xi, cosθ_nodes[idx])
            total += calculate_log_term(E_i, mu_i, T_fm, Φ, Φ̄) * coefficients[idx]
        end
    end
    return -2 * T_fm * total
end

function calculate_pressure(x_state::SVector{5, TF}, mu_vec::AbstractVector{TM}, T_fm::TR, thermal_nodes, xi) where {TF, TM, TR}
    φ = SVector{3, TF}(x_state[1], x_state[2], x_state[3])
    Φ, Φ̄ = x_state[4], x_state[5]

    chi = calculate_chiral(φ)
    U = calculate_U(T_fm, Φ, Φ̄)
    masses = calculate_mass_vec(φ)

    thermal_p_mesh, cosθ_mesh, thermal_coefficients = thermal_nodes

    energy_sum = calculate_energy_sum(masses)
    log_sum = calculate_log_sum(masses, thermal_p_mesh, cosθ_mesh, thermal_coefficients, Φ, Φ̄, mu_vec, T_fm, xi)
    return -(chi + U + energy_sum + log_sum)
end

"""残差方程组的核心，计算偏导数向量，解方程意味着这些偏导数为零"""
function calculate_core(x_state::SVector{5, TF}, mu_vec::AbstractVector{TM}, T_fm::TR, thermal_nodes, xi) where {TF, TM, TR}
    # keep mu_vec as-is (may be Float64 or Dual). Numeric operations will promote as needed
    pressure_fn = y -> begin
        eltp = typeof(y[1])
        y_s = SVector{5, eltp}(Tuple(y))
        calculate_pressure(y_s, mu_vec, T_fm, thermal_nodes, xi)
    end
    grad = ForwardDiff.gradient(pressure_fn, x_state)
    grad_type = typeof(grad[1])
    return SVector{5, grad_type}(Tuple(grad))
end

"""计算粒子数密度ρ"""
function calculate_rho(x_state::SVector{5, TF}, mu_vec::AbstractVector{TM}, T_fm::TR, thermal_nodes, xi) where {TF, TM, TR}
    # mu_vec may be Float64 or Dual; let ForwardDiff handle promotion when computing gradient
    pressure_mu = μ -> calculate_pressure(x_state, μ, T_fm, thermal_nodes, xi)
    grad = ForwardDiff.gradient(pressure_mu, mu_vec)
    grad_type = typeof(grad[1])
    return SVector{3, grad_type}(Tuple(grad))
end

"""计算所有热力学量，返回压力、归一化密度、熵、能量"""
function calculate_thermo(x_state::SVector{5, TF}, mu_vec::AbstractVector{TM}, T_fm::TR, thermal_nodes, xi) where {TF, TM, TR}
    # keep mu_vec in its incoming element type; calculus routines will promote as needed
    rho_vec = calculate_rho(x_state, mu_vec, T_fm, thermal_nodes, xi)
    rho_norm = sum(rho_vec) / (3.0 * ρ0)
    pressure_T = τ -> calculate_pressure(x_state, mu_vec, τ, thermal_nodes, xi)
    entropy = ForwardDiff.derivative(pressure_T, T_fm)
    pressure = calculate_pressure(x_state, mu_vec, T_fm, thermal_nodes, xi)
    # ensure mu_vec elements are compatible for multiplication with rho_vec (promotion will occur)
    energy = -pressure + sum(mu_vec .* rho_vec) + T_fm * entropy
    return pressure, rho_norm, entropy, energy
end

"""T-rho确定时的残差方程组"""
function residual_rho!(fvec, x, T_fm, rho_target, thermal_nodes, xi)
    eltp = typeof(x[1])
    x_state = SVector{5, eltp}(Tuple(x[1:5]))
    mu_state = SVector{3, eltp}(x[6], x[7], x[8])
    fvec[1:5] = calculate_core(x_state, mu_state, T_fm, thermal_nodes, xi)
    fvec[6] = x[6] - x[7]
    fvec[7] = x[7] - x[8]
    rho = calculate_rho(x_state, mu_state, T_fm, thermal_nodes, xi)
    fvec[8] = sum(rho) / (3.0 * ρ0) - rho_target
    return nothing
end

"""T-mu确定时的残差方程组"""
function residual_mu!(fvec, x, mu_vec, T_fm, thermal_nodes, xi)
    eltp = typeof(x[1])
    x_state = SVector{5, eltp}(Tuple(x))
    core_grad = calculate_core(x_state, mu_vec, T_fm, thermal_nodes, xi)
    fvec .= core_grad
    return nothing
end

function solve_fixed_rho(T_mev::Float64, rho_target::Float64; xi::Float64=0.0, seed_state::AbstractVector=DEFAULT_RHO_GUESS, p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT, nlsolve_kwargs...)
    T_fm = T_mev / ħc_MeV_fm
    thermal_nodes = cached_nodes(p_num, t_num)
    x0 = length(seed_state) == 8 ? Float64.(seed_state) : Float64.(DEFAULT_RHO_GUESS)
    f! = (F, x) -> residual_rho!(F, x, T_fm, rho_target, thermal_nodes, xi)
    extra = (; nlsolve_kwargs...)
    res = nlsolve(f!, x0; autodiff=:forward, method=:newton, xtol=1e-9, ftol=1e-9, extra...)
    x_sol = res.zero
    x_state = SVector{5}(Tuple(x_sol[1:5]))
    mu_state = SVector{3}(x_sol[6], x_sol[7], x_sol[8])
    pressure, rho_norm, entropy, energy = calculate_thermo(x_state, mu_state, T_fm, thermal_nodes, xi)
    return SolverResult(
        :rho,
        res.f_converged,
        Vector{Float64}(x_sol),
        mu_state,
        pressure,
        rho_norm,
        entropy,
        energy,
        res.iterations,
        res.residual_norm,
        xi,
    )
end

function solve_fixed_mu(T_mev::Float64, mu_fm::Float64; xi::Float64=0.0, seed_state::AbstractVector=DEFAULT_MU_GUESS, p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT, nlsolve_kwargs...)
    T_fm = T_mev / ħc_MeV_fm
    thermal_nodes = cached_nodes(p_num, t_num)
    x0 = length(seed_state) == 5 ? Float64.(seed_state) : Float64.(DEFAULT_MU_GUESS)
    mu_vec = SVector{3}(mu_fm, mu_fm, mu_fm)
    f! = (F, x) -> residual_mu!(F, x, mu_vec, T_fm, thermal_nodes, xi)
    extra = (; nlsolve_kwargs...)
    res = nlsolve(f!, x0; autodiff=:forward, method=:newton, xtol=1e-9, ftol=1e-9, extra...)
    x_state = SVector{5}(Tuple(res.zero))
    pressure, rho_norm, entropy, energy = calculate_thermo(x_state, mu_vec, T_fm, thermal_nodes, xi)
    return SolverResult(
        :mu,
        res.f_converged,
        Vector{Float64}(res.zero),
        mu_vec,
        pressure,
        rho_norm,
        entropy,
        energy,
        res.iterations,
        res.residual_norm,
        xi,
    )
end

end # module