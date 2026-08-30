"""
    MagneticIntegrals

外磁场 PNJL 的磁场积分工具。

默认生产路线使用 MFIR：零场三动量截断真空项由 `PNJLCore` 组装，
本模块提供有限磁场的 Hurwitz-zeta 修正；热介质项仍使用 Landau 能级求和。
完整 Landau 真空项保留为显式 legacy/diagnostic 路线，不作为默认正则化。

单位约定：
- 动量/质量/温度/化学势：fm⁻¹
- eB：fm⁻²（即电荷单位 e 已并入）
"""
module MagneticIntegrals

using Base.MathConstants: π
using StaticArrays

const _GAUSSLEGENDRE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "integration", "GaussLegendre.jl"))
if !isdefined(Main, :GaussLegendre)
    Base.include(Main, _GAUSSLEGENDRE_PATH)
end
using Main.GaussLegendre: gauleg

const _CONSTANTS_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "constants", "Constants_PNJL.jl"))
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PATH)
end
using Main.Constants_PNJL: Λ_inv_fm, ħc_MeV_fm, N_color

export QUARK_CHARGE_ABS
export alpha_n, energy_landau, smooth_cutoff
export pz_nodes, resolve_nmax_from_cutoff
export DEFAULT_ZETA_COUNT, zeta_nodes, zeta_prime_minus_one, omega_magnetic_mfir
export omega0_flavor_landau, omegat_flavor_landau
export density_flavor_landau

const QUARK_CHARGE_ABS = SVector{3, Float64}(2 / 3, 1 / 3, 1 / 3)
const MAGNETIC_EB_MIN_MEV2 = 100.0
const MAGNETIC_EB_MIN_FM2 = MAGNETIC_EB_MIN_MEV2 / ħc_MeV_fm^2
const _PZ_NODE_CACHE = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
const DEFAULT_ZETA_COUNT = 64
const _ZETA_NODE_CACHE = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
const _LOG_EPS = 1e-16
const _EXP_LIMIT = 745.0

@inline alpha_n(n::Integer) = n == 0 ? 1.0 : 2.0

"""Validate the positive magnetic-field contract in internal `fm^-2` units."""
@inline function validate_magnetic_eB(eB_fm2::Real)
    value = Float64(eB_fm2)
    isfinite(value) || throw(ArgumentError("magnetic eB_fm2 must be finite, got $(eB_fm2)"))
    value >= MAGNETIC_EB_MIN_FM2 || throw(ArgumentError(
        "magnetic eB_fm2 must be >= $(MAGNETIC_EB_MIN_FM2) " *
        "(equivalent to eB >= $(MAGNETIC_EB_MIN_MEV2) MeV^2), got $(eB_fm2)",
    ))
    return value
end

@inline function _energy_landau_unchecked(mass::Real, pz::Real, n::Integer, q_abs::Real, eB::Real)
    return sqrt(2 * n * q_abs * eB + pz^2 + mass^2)
end

@inline function energy_landau(mass::Real, pz::Real, n::Integer, q_abs::Real, eB::Real)
    return _energy_landau_unchecked(mass, pz, n, q_abs, validate_magnetic_eB(eB))
end

@inline function smooth_cutoff(p::Real; Λ::Real=Λ_inv_fm, N::Int=10)
    N >= 1 || throw(ArgumentError("N must be >= 1, got $N"))
    ratio = (abs(p) / max(abs(Λ), 1e-12))^(2N)
    return inv(sqrt(1 + ratio))
end

function pz_nodes(p_num::Int)
    p_num >= 4 || throw(ArgumentError("p_num must be >= 4, got $p_num"))
    return get!(_PZ_NODE_CACHE, p_num) do
        nodes, weights = gauleg(0.0, 1.0, p_num)
        (nodes, weights)
    end
end

function zeta_nodes(zeta_num::Int=DEFAULT_ZETA_COUNT)
    zeta_num >= 8 || throw(ArgumentError("zeta_num must be >= 8, got $(zeta_num)"))
    return get!(_ZETA_NODE_CACHE, zeta_num) do
        gauleg(0.0, 1.0, zeta_num)
    end
end

# Jet-compatible Abel-Plana representation of zeta'(-1, x). The fixed
# quadrature nodes are numerical data; operations on x remain generic so
# ForwardDiff can differentiate the MFIR correction with respect to the gap.
@inline function _atan_series_small(x::Real; terms::Int=18)
    x_squared = x * x
    power = x
    total = zero(x)
    sign = one(x)
    @inbounds for index in 0:(terms - 1)
        total += sign * power / (2 * index + 1)
        power *= x_squared
        sign = -sign
    end
    return total
end

@inline function _atan_nonnegative(x::Real)
    reduced = x
    @inbounds for _ in 1:2
        reduced /= one(reduced) + sqrt(one(reduced) + reduced * reduced)
    end
    return 4 * _atan_series_small(reduced)
end

function zeta_prime_minus_one(x::Real; zeta_num::Int=DEFAULT_ZETA_COUNT)
    isfinite(x) && x > zero(x) || throw(ArgumentError(
        "Hurwitz-zeta argument x must be finite and positive, got $(x)",
    ))
    nodes, weights = zeta_nodes(zeta_num)
    scalar_term = x * (-x + 2 * (x - one(x)) * log(x)) / 4
    integral = zero(x)
    @inbounds for i in eachindex(nodes, weights)
        u = nodes[i]
        momentum = -log1p(-u)
        jacobian = inv(one(u) - u)
        integrand = (
            2 * x * _atan_nonnegative(momentum / x) +
            momentum * log(x * x + momentum * momentum)
        ) / (2 * (exp(2 * π * momentum) - 1))
        integral += weights[i] * integrand * jacobian
    end
    return scalar_term + 2 * integral
end

"""Finite MFIR magnetic-vacuum correction for one flavor.

The zero-field cutoff vacuum term is deliberately not included here. The
caller adds it once using the shared `PNJLCore.vacuum_integral_with_cutoff`,
then adds this finite Hurwitz-zeta correction.
"""
function omega_magnetic_mfir(
    mass::Real,
    q_abs::Real,
    eB::Real;
    zeta_num::Int=DEFAULT_ZETA_COUNT,
    N_c::Real=N_color,
)
    eB_value = validate_magnetic_eB(eB)
    qB = q_abs * eB_value
    mass_abs = abs(mass)
    mass_safe = mass_abs + one(mass_abs) * 1e-12
    x = mass_safe * mass_safe / (2 * qB)
    bracket = zeta_prime_minus_one(x; zeta_num=zeta_num) -
        (x * x - x) * log(x) / 2 + x * x / 4
    return -N_c * qB * qB * bracket / (2 * π^2)
end

@inline function _safe_log(x)
    return log(max(float(x), _LOG_EPS))
end

@inline function _scaled_polyakov_terms(x::Real, Φ::Real, Φbar::Real)
    # Scale the four polynomial terms before exponentiating. This avoids
    # overflow at low T while preserving both the log and its mu derivative.
    m = max(0.0, x, 2 * x, 3 * x)
    e1 = exp(clamp(x - m, -_EXP_LIMIT, 0.0))
    e2 = exp(clamp(2 * x - m, -_EXP_LIMIT, 0.0))
    e3 = exp(clamp(3 * x - m, -_EXP_LIMIT, 0.0))
    terms = (exp(-m), 3 * Φ * e1, 3 * Φbar * e2, e3)
    return m, terms
end

@inline function _polyakov_log_and_net_density(E::Real, mu::Real, T::Real, Φ::Real, Φbar::Real)
    a = -(E - mu) / T
    b = -(E + mu) / T
    ma, plus = _scaled_polyakov_terms(a, Φ, Φbar)
    mb, minus = _scaled_polyakov_terms(b, Φbar, Φ)
    f_plus = max(sum(plus), _LOG_EPS)
    f_minus = max(sum(minus), _LOG_EPS)
    plus_mu = (plus[2] + 2 * plus[3] + 3 * plus[4]) / f_plus
    minus_mu = (minus[2] + 2 * minus[3] + 3 * minus[4]) / f_minus
    return ma + log(f_plus) + mb + log(f_minus), plus_mu - minus_mu
end

@inline function _log_polyakov_pair(E::Real, mu::Real, T::Real, Φ::Real, Φbar::Real)
    return first(_polyakov_log_and_net_density(E, mu, T, Φ, Φbar))
end

function resolve_nmax_from_cutoff(mass::Real, mu::Real, q_abs::Real, eB::Real; Λ::Real=Λ_inv_fm)
    eB_value = validate_magnetic_eB(eB)
    p2_eff = max(Λ^2, mu^2) + mass^2
    nmax = floor(Int, p2_eff / (2 * q_abs * eB_value))
    return max(nmax, 0)
end

@inline function _integrate_pz_even(f, p_num::Int, pz_max::Real)
    nodes, weights = pz_nodes(p_num)
    scale = pz_max
    acc = zero(Float64)
    @inbounds for i in eachindex(nodes)
        pz = scale * nodes[i]
        acc += weights[i] * f(pz)
    end
    return 2 * scale * acc
end

function omega0_flavor_landau(
    mass::Real,
    q_abs::Real,
    eB::Real;
    n_max::Int,
    p_num::Int=96,
    pz_max::Real=max(8 * Λ_inv_fm, 25.0),
    cutoff_N::Int=10,
)
    eB_value = validate_magnetic_eB(eB)
    pref = -3.0 * q_abs * eB_value / (2 * π)
    total = 0.0
    @inbounds for n in 0:n_max
        an = alpha_n(n)
        int_val = _integrate_pz_even(p_num, pz_max) do pz
            E = _energy_landau_unchecked(mass, pz, n, q_abs, eB_value)
            p3 = sqrt(pz^2 + 2 * n * q_abs * eB_value)
            fc = smooth_cutoff(p3; N=cutoff_N)
            fc^2 * E
        end
        total += an * int_val / (2 * π)
    end
    return pref * total
end

function omegat_flavor_landau(
    mass::Real,
    mu::Real,
    T::Real,
    Φ::Real,
    Φbar::Real,
    q_abs::Real,
    eB::Real;
    n_max::Int,
    p_num::Int=96,
    pz_max::Real=max(8 * Λ_inv_fm, 25.0),
)
    eB_value = validate_magnetic_eB(eB)
    T <= 1e-12 && return 0.0

    pref = -T * q_abs * eB_value / (2 * π)
    total = 0.0
    @inbounds for n in 0:n_max
        an = alpha_n(n)
        int_val = _integrate_pz_even(p_num, pz_max) do pz
            E = _energy_landau_unchecked(mass, pz, n, q_abs, eB_value)
            _log_polyakov_pair(E, mu, T, Φ, Φbar)
        end
        total += an * int_val / (2 * π)
    end
    return pref * total
end

function density_flavor_landau(
    mass::Real,
    mu::Real,
    T::Real,
    Φ::Real,
    Φbar::Real,
    q_abs::Real,
    eB::Real;
    n_max::Int,
    p_num::Int=96,
    pz_max::Real=max(8 * Λ_inv_fm, 25.0),
)
    eB_value = validate_magnetic_eB(eB)
    T <= 1e-12 && return 0.0

    # The Polyakov polynomial already contains the color trace. Its mu
    # derivative therefore carries the color factor and the prefactor is the
    # Landau phase-space measure only.
    pref = q_abs * eB_value / (2 * π)
    total = 0.0
    @inbounds for n in 0:n_max
        an = alpha_n(n)
        int_val = _integrate_pz_even(p_num, pz_max) do pz
            E = _energy_landau_unchecked(mass, pz, n, q_abs, eB_value)
            _, net = _polyakov_log_and_net_density(E, mu, T, Φ, Φbar)
            net
        end
        total += an * int_val / (2 * π)
    end
    return pref * total
end

end # module MagneticIntegrals
