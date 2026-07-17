"""PNJLIntegrals

PNJL 模型在 models 路径下使用的最小热项积分与节点缓存实现。

目标：让 `models/pnjl_physics/PNJLModel.jl` 的热项计算不再依赖 legacy 的
`src/pnjl/core/Integrals.jl`（以及其上层 Thermodynamics include）。

说明：此处实现保持与 legacy 积分节点与 log-sum 形式一致，但刻意只提供
models 侧当前所需的最小 API：`cached_nodes` 与 `calculate_log_sum`。
"""
module PNJLIntegrals

using Base.MathConstants: π
using ForwardDiff
using QuadGK: quadgk
using StaticArrays

# Gauss-Legendre nodes/weights
const _GAUSSLEGENDRE_PATH = normpath(joinpath(@__DIR__, "..", "..", "integration", "GaussLegendre.jl"))
if !isdefined(Main, :GaussLegendre)
    Base.include(Main, _GAUSSLEGENDRE_PATH)
end

using Main.GaussLegendre:
    gauleg,
    DEFAULT_COSΘ_HALF_NODES,
    DEFAULT_COSΘ_HALF_WEIGHTS,
    DEFAULT_MOMENTUM_NODES,
    DEFAULT_MOMENTUM_WEIGHTS

export cached_nodes, calculate_log_sum
export calculate_log_sum_rs_reduced_adaptive
export calculate_log_sum_rs_reduced_adaptive_with_error
export integrate_rs_reduced_radial
export rs_anisotropy_measure_factor
export validate_rs_anisotropy
export validate_thermal_quadrature_policy
export validate_thermal_quadrature_controls
export SUPPORTED_THERMAL_QUADRATURE_POLICIES
export DEFAULT_THETA_COUNT, DEFAULT_MOMENTUM_COUNT, DEFAULT_THERMAL_P_MAX_INV_FM

const DEFAULT_THETA_COUNT = length(DEFAULT_COSΘ_HALF_NODES)
const DEFAULT_MOMENTUM_COUNT = length(DEFAULT_MOMENTUM_NODES)
const DEFAULT_THERMAL_P_MAX_INV_FM = 10.0
const SUPPORTED_THERMAL_QUADRATURE_POLICIES = (:tensor_gauss, :rs_reduced_adaptive)
const THETA_DEFAULT_NODES = DEFAULT_COSΘ_HALF_NODES
const THETA_DEFAULT_WEIGHTS = DEFAULT_COSΘ_HALF_WEIGHTS .* 2.0
const THERMAL_DEFAULT_NODES = DEFAULT_MOMENTUM_NODES
const THERMAL_DEFAULT_WEIGHTS = DEFAULT_MOMENTUM_WEIGHTS

"""积分节点缓存：(p_num, t_num, p_max_inv_fm) -> (p_mesh, cosθ_mesh, coefficients)"""
const NODE_CACHE = Dict{Tuple{Int, Int, Float64}, NTuple{3, Matrix{Float64}}}()

@inline function validate_thermal_quadrature_policy(policy::Symbol)
    policy in SUPPORTED_THERMAL_QUADRATURE_POLICIES || throw(ArgumentError(
        "unsupported thermo_quadrature_policy=$(policy); expected one of $(SUPPORTED_THERMAL_QUADRATURE_POLICIES)",
    ))
    return policy
end

@inline function validate_thermal_quadrature_controls(
    rtol::Float64,
    atol::Float64,
    maxevals::Int,
)
    isfinite(rtol) && rtol > 0 || throw(ArgumentError("thermo quadrature rtol must be finite and positive, got $(rtol)"))
    isfinite(atol) && atol >= 0 || throw(ArgumentError("thermo quadrature atol must be finite and nonnegative, got $(atol)"))
    maxevals > 0 || throw(ArgumentError("thermo quadrature maxevals must be positive, got $(maxevals)"))
    return nothing
end

@inline _primal_float(x) = Float64(x)
@inline _primal_float(x::ForwardDiff.Dual) = _primal_float(ForwardDiff.value(x))

@inline function validate_rs_anisotropy(xi)
    xi_value = _primal_float(xi)
    isfinite(xi_value) || throw(ArgumentError("xi must be finite, got $(xi_value)"))
    xi_value > -1.0 || throw(ArgumentError("RS thermodynamic angular reduction requires xi > -1, got $(xi_value)"))
    return xi
end

@inline function rs_anisotropy_measure_factor(xi)
    validate_rs_anisotropy(xi)
    return inv(sqrt(one(xi) + xi))
end

@inline function theta_nodes(t_num::Int)
    if t_num == DEFAULT_THETA_COUNT
        return THETA_DEFAULT_NODES, THETA_DEFAULT_WEIGHTS
    end
    nodes, weights = gauleg(0.0, 1.0, t_num)
    return nodes, weights .* 2.0
end

@inline function thermal_nodes(p_num::Int, p_max_inv_fm::Float64=DEFAULT_THERMAL_P_MAX_INV_FM)
    if p_num == DEFAULT_MOMENTUM_COUNT && isapprox(p_max_inv_fm, DEFAULT_THERMAL_P_MAX_INV_FM; atol=1e-12, rtol=0.0)
        return THERMAL_DEFAULT_NODES, THERMAL_DEFAULT_WEIGHTS
    end
    return gauleg(0.0, p_max_inv_fm, p_num)
end

function build_nodes(p_num::Int, t_num::Int, p_max_inv_fm::Float64=DEFAULT_THERMAL_P_MAX_INV_FM)
    momentum_nodes, momentum_weights = thermal_nodes(p_num, p_max_inv_fm)
    cosθ_nodes, cosθ_weights = theta_nodes(t_num)

    thermal_p_mesh = repeat(momentum_nodes, 1, t_num)
    cosθ_mesh = repeat(cosθ_nodes', p_num, 1)
    weight_mesh = momentum_weights * cosθ_weights'
    thermal_coefficients = weight_mesh .* thermal_p_mesh .^ 2 ./ (2 * π) ^ 2

    return (thermal_p_mesh, cosθ_mesh, thermal_coefficients)
end

function cached_nodes(p_num::Int, t_num::Int; p_max_inv_fm::Float64=DEFAULT_THERMAL_P_MAX_INV_FM)
    key = (p_num, t_num, p_max_inv_fm)
    return get!(NODE_CACHE, key) do
        build_nodes(p_num, t_num, p_max_inv_fm)
    end
end

const POLYAKOV_EPS = 1e-16

@inline function _safe_log(x; min_val=POLYAKOV_EPS)
    x <= 0 && return log(min_val)
    return x < min_val ? log(min_val) : log(x)
end

@inline function _calculate_rs_distribution_energy(mass_i, p, xi, cosθ)
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

@inline function _fermi_momentum_primal(mass, mu)
    mass_value = abs(_primal_float(mass))
    mu_value = abs(_primal_float(mu))
    return mu_value > mass_value ? sqrt(max(mu_value * mu_value - mass_value * mass_value, 0.0)) : 0.0
end

"""Integrate a radial RS-reduced kernel on `[0, Inf)` with a Fermi-surface breakpoint.

The callback must already include the radial `q^2` measure.  The returned tuple is
`(value, error)` from `QuadGK`.  Breakpoint placement uses primal values only so
ForwardDiff differentiates the physical integrand rather than the numerical node
selection rule.
"""
function integrate_rs_reduced_radial(
    f,
    mass,
    mu;
    rtol::Float64=1e-8,
    atol::Float64=1e-10,
    maxevals::Int=10^7,
)
    validate_thermal_quadrature_controls(rtol, atol, maxevals)

    q_fermi = _fermi_momentum_primal(mass, mu)
    if q_fermi > 0.0
        return quadgk(f, 0.0, q_fermi, Inf; rtol=rtol, atol=atol, maxevals=maxevals)
    end
    return quadgk(f, 0.0, Inf; rtol=rtol, atol=atol, maxevals=maxevals)
end

@inline function _zero_temperature_pressure_radial(mass, mu_abs)
    mass_abs = abs(mass)
    mass_value = abs(_primal_float(mass_abs))
    mu_value = _primal_float(mu_abs)
    mu_value > mass_value || return zero(promote_type(typeof(mass), typeof(mu_abs)))

    q_fermi = sqrt(mu_abs * mu_abs - mass_abs * mass_abs)
    if mass_value <= 1e-14
        return mu_abs^4 / 12
    end

    energy_integral = (
        q_fermi * mu_abs * (2 * q_fermi * q_fermi + mass_abs * mass_abs) -
        mass_abs^4 * asinh(q_fermi / mass_abs)
    ) / 8
    return mu_abs * q_fermi^3 / 3 - energy_integral
end

function _calculate_log_sum_rs_reduced_zero_temperature(masses, mu_vec, xi)
    TT = promote_type(eltype(masses), eltype(mu_vec), typeof(xi))
    total = zero(TT)
    @inbounds for i in 1:3
        total += _zero_temperature_pressure_radial(
            convert(TT, masses[i]),
            abs(convert(TT, mu_vec[i])),
        )
    end
    radial_measure = convert(TT, 2) * rs_anisotropy_measure_factor(convert(TT, xi)) /
        convert(TT, (2 * π)^2)
    return -convert(TT, 6) * radial_measure * total
end

"""RS-reduced adaptive PNJL thermal grand-potential contribution.

This path applies only to scalar thermodynamic kernels whose angular dependence is
entirely contained in `E_xi`.  It integrates the transformed radial momentum to
infinity, avoiding the finite-spherical-cutoff mismatch of a naive angular
prefactor replacement.
"""
function calculate_log_sum_rs_reduced_adaptive_with_error(
    masses::SVector{3, TF},
    Φ,
    Φ̄,
    mu_vec,
    T_fm,
    xi;
    rtol::Float64=1e-8,
    atol::Float64=1e-10,
    maxevals::Int=10^7,
) where {TF}
    validate_thermal_quadrature_controls(rtol, atol, maxevals)
    T_value = _primal_float(T_fm)
    isfinite(T_value) || throw(ArgumentError("T_fm must be finite, got $(T_value)"))
    T_value >= 0.0 || throw(ArgumentError("T_fm must be nonnegative, got $(T_value)"))
    if T_value == 0.0
        value = _calculate_log_sum_rs_reduced_zero_temperature(masses, mu_vec, xi)
        return (value=value, error=0.0)
    end

    TT = promote_type(TF, typeof(Φ), typeof(Φ̄), typeof(T_fm), eltype(mu_vec), typeof(xi))
    ΦT = convert(TT, Φ)
    Φ̄T = convert(TT, Φ̄)
    Tt = convert(TT, T_fm)
    xit = convert(TT, xi)
    radial_measure = convert(TT, 2) * rs_anisotropy_measure_factor(xit) /
        convert(TT, (2 * π)^2)

    total = zero(TT)
    total_error = 0.0
    @inbounds for i in 1:3
        mass_i = convert(TT, masses[i])
        mu_i = convert(TT, mu_vec[i])
        integrand = q -> begin
            qT = convert(TT, q)
            E_dist = sqrt(qT * qT + mass_i * mass_i)
            return qT * qT * (-convert(TT, 2) * Tt) *
                calculate_log_term(E_dist, mu_i, Tt, ΦT, Φ̄T)
        end
        value, error = integrate_rs_reduced_radial(
            integrand,
            mass_i,
            mu_i;
            rtol=rtol,
            atol=atol,
            maxevals=maxevals,
        )
        total += value
        total_error += _primal_float(error)
    end
    return (
        value=radial_measure * total,
        error=abs(_primal_float(radial_measure)) * total_error,
    )
end

function calculate_log_sum_rs_reduced_adaptive(
    masses::SVector{3, TF},
    Φ,
    Φ̄,
    mu_vec,
    T_fm,
    xi;
    rtol::Float64=1e-8,
    atol::Float64=1e-10,
    maxevals::Int=10^7,
) where {TF}
    result = calculate_log_sum_rs_reduced_adaptive_with_error(
        masses,
        Φ,
        Φ̄,
        mu_vec,
        T_fm,
        xi;
        rtol=rtol,
        atol=atol,
        maxevals=maxevals,
    )
    return result.value
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
            E_dist_i = _calculate_rs_distribution_energy(mass_i, p, xit, cosθ)
            total += calculate_log_term(E_dist_i, mu_i, Tt, ΦT, Φ̄T) * w
        end
    end

    return -convert(TT, 2) * Tt * total
end

end # module PNJLIntegrals
