"""
    MagneticIntegrals

外磁场 PNJL 的 Landau 能级积分工具。

单位约定：
- 动量/质量/温度/化学势：fm⁻¹
- eB：fm⁻²（即电荷单位 e 已并入）
"""
module MagneticIntegrals

using Base.MathConstants: π
using StaticArrays

const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

const _GAUSSLEGENDRE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "integration", "GaussLegendre.jl"))
IncludeOnce.include_once!(Main, :GaussLegendre, _GAUSSLEGENDRE_PATH)
using Main.GaussLegendre: gauleg

const _CONSTANTS_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "constants", "Constants_PNJL.jl"))
IncludeOnce.include_once!(Main, :Constants_PNJL, _CONSTANTS_PATH)
using Main.Constants_PNJL: Λ_inv_fm

export QUARK_CHARGE_ABS
export alpha_n, energy_landau, smooth_cutoff
export pz_nodes, resolve_nmax_from_cutoff
export omega0_flavor_landau, omegat_flavor_landau
export density_flavor_landau

const QUARK_CHARGE_ABS = SVector{3, Float64}(2 / 3, 1 / 3, 1 / 3)
const _PZ_NODE_CACHE = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
const _LOG_EPS = 1e-16

@inline alpha_n(n::Integer) = n == 0 ? 1.0 : 2.0

@inline function energy_landau(mass::Real, pz::Real, n::Integer, q_abs::Real, eB::Real)
    return sqrt(2 * n * abs(q_abs * eB) + pz^2 + mass^2)
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

@inline function _safe_log(x)
    return log(max(float(x), _LOG_EPS))
end

@inline function _log_polyakov_pair(E::Real, mu::Real, T::Real, Φ::Real, Φbar::Real)
    a = -(E - mu) / T
    b = -(E + mu) / T
    e1a = exp(a)
    e2a = e1a * e1a
    e3a = e2a * e1a
    e1b = exp(b)
    e2b = e1b * e1b
    e3b = e2b * e1b
    f_plus = 1 + 3 * Φ * e1a + 3 * Φbar * e2a + e3a
    f_minus = 1 + 3 * Φbar * e1b + 3 * Φ * e2b + e3b
    return _safe_log(f_plus) + _safe_log(f_minus)
end

function resolve_nmax_from_cutoff(mass::Real, mu::Real, q_abs::Real, eB::Real; Λ::Real=Λ_inv_fm)
    abs(q_abs * eB) <= 1e-14 && return 0
    p2_eff = max(Λ^2, mu^2) + mass^2
    nmax = floor(Int, p2_eff / (2 * abs(q_abs * eB)))
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
    abs(q_abs * eB) <= 1e-14 && return 0.0
    pref = -3.0 * abs(q_abs * eB) / (2 * π)
    total = 0.0
    @inbounds for n in 0:n_max
        an = alpha_n(n)
        int_val = _integrate_pz_even(p_num, pz_max) do pz
            E = energy_landau(mass, pz, n, q_abs, eB)
            p3 = sqrt(pz^2 + 2 * n * abs(q_abs * eB))
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
    abs(q_abs * eB) <= 1e-14 && return 0.0
    T <= 1e-12 && return 0.0

    pref = -T * abs(q_abs * eB) / (2 * π)
    total = 0.0
    @inbounds for n in 0:n_max
        an = alpha_n(n)
        int_val = _integrate_pz_even(p_num, pz_max) do pz
            E = energy_landau(mass, pz, n, q_abs, eB)
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
    T <= 1e-12 && return 0.0
    abs(q_abs * eB) <= 1e-14 && return 0.0

    pref = 3.0 * abs(q_abs * eB) / (2 * π)
    total = 0.0
    @inbounds for n in 0:n_max
        an = alpha_n(n)
        int_val = _integrate_pz_even(p_num, pz_max) do pz
            E = energy_landau(mass, pz, n, q_abs, eB)
            x = (E - mu) / T
            y = (E + mu) / T
            fq = inv(1 + exp(3 * x))
            fa = inv(1 + exp(3 * y))
            fq - fa
        end
        total += an * int_val / (2 * π)
    end
    return pref * total
end

end # module MagneticIntegrals

