"""PNJLIntegrals

PNJL 模型在 models 路径下使用的最小热项积分与节点缓存实现。

目标：让 `src/models/pnjl/PNJLModel.jl` 的热项计算不再依赖 legacy 的
`src/pnjl/core/Integrals.jl`（以及其上层 Thermodynamics include）。

说明：此处实现保持与 legacy 积分节点与 log-sum 形式一致，但刻意只提供
models 侧当前所需的最小 API：`cached_nodes` 与 `calculate_log_sum`。
"""
module PNJLIntegrals

using Base.MathConstants: π
using StaticArrays

# Include-once helper
const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

# Gauss-Legendre nodes/weights
const _GAUSSLEGENDRE_PATH = normpath(joinpath(@__DIR__, "..", "..", "integration", "GaussLegendre.jl"))
IncludeOnce.include_once!(Main, :GaussLegendre, _GAUSSLEGENDRE_PATH)

using Main.GaussLegendre:
    gauleg,
    DEFAULT_COSΘ_HALF_NODES,
    DEFAULT_COSΘ_HALF_WEIGHTS,
    DEFAULT_MOMENTUM_NODES,
    DEFAULT_MOMENTUM_WEIGHTS

export cached_nodes, calculate_log_sum
export DEFAULT_THETA_COUNT, DEFAULT_MOMENTUM_COUNT

const DEFAULT_THETA_COUNT = length(DEFAULT_COSΘ_HALF_NODES)
const DEFAULT_MOMENTUM_COUNT = length(DEFAULT_MOMENTUM_NODES)
const THETA_DEFAULT_NODES = DEFAULT_COSΘ_HALF_NODES
const THETA_DEFAULT_WEIGHTS = DEFAULT_COSΘ_HALF_WEIGHTS .* 2.0
const THERMAL_DEFAULT_NODES = DEFAULT_MOMENTUM_NODES
const THERMAL_DEFAULT_WEIGHTS = DEFAULT_MOMENTUM_WEIGHTS

"""积分节点缓存：(p_num, t_num) -> (p_mesh, cosθ_mesh, coefficients)"""
const NODE_CACHE = Dict{Tuple{Int, Int}, NTuple{3, Matrix{Float64}}}()

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
    thermal_coefficients = weight_mesh .* thermal_p_mesh .^ 2 ./ (2 * π) ^ 2

    return (thermal_p_mesh, cosθ_mesh, thermal_coefficients)
end

function cached_nodes(p_num::Int, t_num::Int)
    key = (p_num, t_num)
    return get!(NODE_CACHE, key) do
        build_nodes(p_num, t_num)
    end
end

const POLYAKOV_EPS = 1e-16

@inline function _safe_log(x; min_val=POLYAKOV_EPS)
    x <= 0 && return log(min_val)
    return x < min_val ? log(min_val) : log(x)
end

@inline function _calculate_energy_anisotropic(mass_i, p, xi, cosθ)
    return sqrt(p ^ 2 + mass_i ^ 2 + xi * (p * cosθ) ^ 2)
end

@inline function calculate_log_term(E_i, mu_i, T_fm, Φ, Φ̄)
    TT = promote_type(typeof(E_i), typeof(mu_i), typeof(T_fm), typeof(Φ), typeof(Φ̄))
    E = convert(TT, E_i)
    μ = convert(TT, mu_i)
    T = convert(TT, T_fm)
    ΦT = convert(TT, Φ)
    Φ̄T = convert(TT, Φ̄)

    invT = one(TT) / T
    epsT = convert(TT, POLYAKOV_EPS)
    three = convert(TT, 3)
    two = convert(TT, 2)

    a = -(E - μ) * invT
    b = -(E + μ) * invT

    if a > 0
        m_a = three * a
        exp_a_m = exp(-two * a)
        exp_2a_m = exp(-a)
        exp_neg_m = exp(-m_a)
        term_a = exp_neg_m + three * ΦT * exp_a_m + three * Φ̄T * exp_2a_m + one(TT)
        log_f_plus = m_a + log(max(term_a, epsT))
    else
        exp_a = exp(a)
        exp_2a = exp_a * exp_a
        exp_3a = exp_a * exp_2a
        f_plus = one(TT) + three * ΦT * exp_a + three * Φ̄T * exp_2a + exp_3a
        log_f_plus = log(max(f_plus, epsT))
    end

    if b > 0
        m_b = three * b
        exp_b_m = exp(-two * b)
        exp_2b_m = exp(-b)
        exp_neg_m = exp(-m_b)
        term_b = exp_neg_m + three * Φ̄T * exp_b_m + three * ΦT * exp_2b_m + one(TT)
        log_f_minus = m_b + log(max(term_b, epsT))
    else
        exp_b = exp(b)
        exp_2b = exp_b * exp_b
        exp_3b = exp_b * exp_2b
        f_minus = one(TT) + three * Φ̄T * exp_b + three * ΦT * exp_2b + exp_3b
        log_f_minus = log(max(f_minus, epsT))
    end

    return log_f_plus + log_f_minus
end

"""calculate_log_sum(masses, p_mesh, cosθ_mesh, coefficients, Φ, Φ̄, mu_vec, T_fm, xi) -> Float64

热项对数和（各向异性 RS 形式），保持与 legacy 版本一致：

-2T ∑_i ∫ d^3p/(2π)^3 [Q1 + Q2]

注意：当前实现以 Float64 为主，供 models 侧过渡期使用。
"""
function calculate_log_sum(
    masses::SVector{3, TF},
    p_mesh::AbstractMatrix,
    cosθ_mesh::AbstractMatrix,
    coefficients::AbstractMatrix,
    Φ,
    Φ̄,
    mu_vec,
    T_fm,
    xi,
) where {TF}
    TT = promote_type(TF, typeof(Φ), typeof(Φ̄), typeof(T_fm), eltype(mu_vec), typeof(xi))

    ΦT = convert(TT, Φ)
    Φ̄T = convert(TT, Φ̄)
    Tt = convert(TT, T_fm)
    xit = convert(TT, xi)

    total = zero(TT)

    @inbounds for i in 1:3
        mass_i = convert(TT, masses[i])
        mu_i = convert(TT, mu_vec[i])
        for idx in eachindex(p_mesh)
            p = convert(TT, p_mesh[idx])
            cosθ = convert(TT, cosθ_mesh[idx])
            w = convert(TT, coefficients[idx])
            E_i = _calculate_energy_anisotropic(mass_i, p, xit, cosθ)
            total += calculate_log_term(E_i, mu_i, Tt, ΦT, Φ̄T) * w
        end
    end

    return -convert(TT, 2) * Tt * total
end

end # module PNJLIntegrals
