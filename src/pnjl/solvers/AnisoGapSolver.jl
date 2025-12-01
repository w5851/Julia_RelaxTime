module AnisoGapSolver

using LinearAlgebra
using Base.MathConstants: π
using StaticArrays
using ForwardDiff
using NLsolve
using SpecialFunctions: log
using FastGaussQuadrature: gausslegendre

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

const ρ0 = ρ0_inv_fm3
const DEFAULT_RHO_GUESS = [-0.07441, -0.07441, -1.86717, 0.02032, 0.02372, 2.15664, 2.15664, 2.15664]
const DEFAULT_MU_GUESS = [-1.8, -1.8, -2.1, 0.8, 0.8]
const NODE_CACHE = Dict{Tuple{Int, Int}, Tuple{NTuple{3, Matrix{Float64}}, NTuple{3, Matrix{Float64}}}}()

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

function gauleg(a::Float64, b::Float64, n::Int)
    nodes, weights = gausslegendre(n)
    scale = (b - a) / 2
    shift = (a + b) / 2
    nodes = @. scale * nodes + shift
    weights = @. scale * weights
    return nodes, weights
end

function build_nodes(p_num::Int, t_num::Int)
    p_nodes1, p_weights1 = gauleg(0.0, Λ_inv_fm, p_num)
    p_nodes2, p_weights2 = gauleg(0.0, 20.0, p_num)
    t_nodes, t_weights = gauleg(0.0, 1.0, t_num)

    p1_mesh = repeat(p_nodes1, 1, t_num)
    t_mesh = repeat(t_nodes', p_num, 1)
    w1_mesh = p_weights1 * t_weights'
    coefficient1 = w1_mesh .* p1_mesh.^2 ./ π^2

    p2_mesh = repeat(p_nodes2, 1, t_num)
    w2_mesh = p_weights2 * t_weights'
    coefficient2 = w2_mesh .* p2_mesh.^2 ./ π^2

    nodes_1 = (p1_mesh, t_mesh, coefficient1)
    nodes_2 = (p2_mesh, t_mesh, coefficient2)
    return nodes_1, nodes_2
end

function cached_nodes(p_num::Int, t_num::Int)
    key = (p_num, t_num)
    return get!(NODE_CACHE, key) do
        build_nodes(p_num, t_num)
    end
end

@inline function calculate_chiral(phi::SVector{3, T}) where {T}
    2 * G_fm2 * sum(phi .^ 2) - 4 * K_fm5 * prod(phi)
end

@inline function calculate_U(T::Real, Phi1::Real, Phi2::Real)
    T_ratio = T0_inv_fm / T
    Ta = a0 + a1 * T_ratio + a2 * T_ratio^2
    Tb = b3 * T_ratio^3
    value = 1 - 6 * Phi2 * Phi1 + 4 * (Phi2^3 + Phi1^3) - 3 * (Phi2 * Phi1)^2
    return T^4 * (-0.5 * Ta * Phi2 * Phi1 + Tb * safe_log(value))
end

@inline function calculate_mass_vec(phi::SVector{3, T}) where {T}
    phiu, phid, phis = phi
    return SVector{3, T}(
        m_ud0_inv_fm - 4 * G_fm2 * phiu + 2 * K_fm5 * phid * phis,
        m_ud0_inv_fm - 4 * G_fm2 * phid + 2 * K_fm5 * phiu * phis,
        m_s0_inv_fm - 4 * G_fm2 * phis + 2 * K_fm5 * phiu * phid,
    )
end

@inline function calculate_energy(mass_i, p, xi, t)
    return sqrt(p^2 + mass_i^2 + xi * (p * t)^2)
end

@inline function calculate_log_term(E_i, mu_i, T, Phi1, Phi2)
    invT = 1.0 / T
    x_i = (E_i - mu_i) * invT
    x_i_anti = (E_i + mu_i) * invT

    exp1 = exp(-x_i)
    exp2 = exp1 * exp1
    exp3 = exp1 * exp2
    exp1_anti = exp(-x_i_anti)
    exp2_anti = exp1_anti * exp1_anti
    exp3_anti = exp1_anti * exp2_anti

    f1_val = 1.0 + 3.0 * Phi1 * exp1 + 3.0 * Phi2 * exp2 + exp3
    f2_val = 1.0 + 3.0 * Phi2 * exp1_anti + 3.0 * Phi1 * exp2_anti + exp3_anti

    return safe_log(f1_val) + safe_log(f2_val)
end

function calculate_energy_sum(masses::SVector{3, T}, p_nodes, t_nodes, coefficient, xi) where {T}
    total = zero(T)
    @inbounds for i in 1:3
        mass_i = masses[i]
        for idx in eachindex(p_nodes)
            E = calculate_energy(mass_i, p_nodes[idx], xi, t_nodes[idx])
            total += E * coefficient[idx]
        end
    end
    return -N_color * total
end

function calculate_log_sum(masses::SVector{3, T}, p_nodes, t_nodes, coefficient, Phi1, Phi2, mu_vec, T, xi) where {T}
    total = zero(T)
    @inbounds for i in 1:3
        mass_i = masses[i]
        mu_i = mu_vec[i]
        for idx in eachindex(p_nodes)
            E_i = calculate_energy(mass_i, p_nodes[idx], xi, t_nodes[idx])
            total += calculate_log_term(E_i, mu_i, T, Phi1, Phi2) * coefficient[idx]
        end
    end
    return -T * total
end

function calculate_pressure(x_state::SVector{5, T}, mu_vec::SVector{3, T}, T_val::T, nodes_1, nodes_2, xi) where {T}
    phi = SVector{3, T}(x_state[1], x_state[2], x_state[3])
    Phi1, Phi2 = x_state[4], x_state[5]

    chi = calculate_chiral(phi)
    U = calculate_U(T_val, Phi1, Phi2)
    masses = calculate_mass_vec(phi)

    p_nodes1, t_nodes1, coef1 = nodes_1
    p_nodes2, t_nodes2, coef2 = nodes_2

    energy_sum = calculate_energy_sum(masses, p_nodes1, t_nodes1, coef1, xi)
    log_sum = calculate_log_sum(masses, p_nodes2, t_nodes2, coef2, Phi1, Phi2, mu_vec, T_val, xi)
    return -(chi + U + energy_sum + log_sum)
end

@inline function pressure_wrapper(x_state::SVector{5, T}, mu_vec::SVector{3, T}, T_val::T, nodes_1, nodes_2, xi) where {T}
    calculate_pressure(x_state, mu_vec, T_val, nodes_1, nodes_2, xi)
end

function calculate_core(x_state::SVector{5, T}, mu_vec::SVector{3, T}, T_val::T, nodes_1, nodes_2, xi) where {T}
    f = y -> pressure_wrapper(y, mu_vec, T_val, nodes_1, nodes_2, xi)
    grad = ForwardDiff.gradient(f, x_state)
    return SVector{5, T}(Tuple(grad))
end

function calculate_rho(x_state::SVector{5, T}, mu_vec::SVector{3, T}, T_val::T, nodes_1, nodes_2, xi) where {T}
    f_mu = μ -> pressure_wrapper(x_state, μ, T_val, nodes_1, nodes_2, xi)
    grad = ForwardDiff.gradient(f_mu, mu_vec)
    return SVector{3, T}(Tuple(grad))
end

function calculate_thermo(x_state::SVector{5, T}, mu_vec::SVector{3, T}, T_val::T, nodes_1, nodes_2, xi) where {T}
    rho_vec = calculate_rho(x_state, mu_vec, T_val, nodes_1, nodes_2, xi)
    rho_norm = sum(rho_vec) / (3.0 * ρ0)
    f_T = τ -> pressure_wrapper(x_state, mu_vec, τ, nodes_1, nodes_2, xi)
    entropy = ForwardDiff.derivative(f_T, T_val)
    pressure = pressure_wrapper(x_state, mu_vec, T_val, nodes_1, nodes_2, xi)
    energy = -pressure + sum(mu_vec .* rho_vec) + T_val * entropy
    return pressure, rho_norm, entropy, energy
end

function residual_rho!(fvec, x, T_val, rho_target, nodes_1, nodes_2, xi)
    x_state = SVector{5}(Tuple(x[1:5]))
    mu_state = SVector{3}(x[6], x[7], x[8])
    fvec[1:5] = calculate_core(x_state, mu_state, T_val, nodes_1, nodes_2, xi)
    fvec[6] = x[6] - x[7]
    fvec[7] = x[7] - x[8]
    rho = calculate_rho(x_state, mu_state, T_val, nodes_1, nodes_2, xi)
    fvec[8] = sum(rho) / (3.0 * ρ0) - rho_target
    return nothing
end

function residual_mu!(fvec, x, mu_vec, T_val, nodes_1, nodes_2, xi)
    x_state = SVector{5}(Tuple(x))
    f = calculate_core(x_state, mu_vec, T_val, nodes_1, nodes_2, xi)
    fvec .= f
    return nothing
end

function solve_fixed_rho(T_fm::Float64, rho_target::Float64; xi::Float64=0.0, seed_state::AbstractVector=DEFAULT_RHO_GUESS, p_num::Int=256, t_num::Int=16, nlsolve_kwargs...)
    nodes_1, nodes_2 = cached_nodes(p_num, t_num)
    x0 = length(seed_state) == 8 ? Float64.(seed_state) : Float64.(DEFAULT_RHO_GUESS)
    f! = (F, x) -> residual_rho!(F, x, T_fm, rho_target, nodes_1, nodes_2, xi)
    extra = (; nlsolve_kwargs...)
    res = nlsolve(f!, x0; autodiff=:forward, method=:newton, xtol=1e-9, ftol=1e-9, extra...)
    x_sol = res.zero
    x_state = SVector{5}(Tuple(x_sol[1:5]))
    mu_state = SVector{3}(x_sol[6], x_sol[7], x_sol[8])
    pressure, rho_norm, entropy, energy = calculate_thermo(x_state, mu_state, T_fm, nodes_1, nodes_2, xi)
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
        norm(res.fvec),
        xi,
    )
end

function solve_fixed_mu(T_fm::Float64, mu_fm::Float64; xi::Float64=0.0, seed_state::AbstractVector=DEFAULT_MU_GUESS, p_num::Int=256, t_num::Int=16, nlsolve_kwargs...)
    nodes_1, nodes_2 = cached_nodes(p_num, t_num)
    x0 = length(seed_state) == 5 ? Float64.(seed_state) : Float64.(DEFAULT_MU_GUESS)
    mu_vec = SVector{3}(mu_fm, mu_fm, mu_fm)
    f! = (F, x) -> residual_mu!(F, x, mu_vec, T_fm, nodes_1, nodes_2, xi)
    extra = (; nlsolve_kwargs...)
    res = nlsolve(f!, x0; autodiff=:forward, method=:newton, xtol=1e-9, ftol=1e-9, extra...)
    x_state = SVector{5}(Tuple(res.zero))
    pressure, rho_norm, entropy, energy = calculate_thermo(x_state, mu_vec, T_fm, nodes_1, nodes_2, xi)
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
        norm(res.fvec),
        xi,
    )
end

end # module
*** End Patch