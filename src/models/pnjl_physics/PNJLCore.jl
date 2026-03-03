"""
    PNJLCore

PNJL 模型可复用 core（参数加载 + 单位换算 + 平均场核心公式）。

说明：
- 当前实现默认支持各向异性工作流（`xi` 在热项积分层 `PNJLIntegrals`/`PNJLDistributions` 处理）。
- 本模块聚焦“参数与局域势能公式”，便于 `PNJLModel` 原生实现，避免依赖 legacy core include 链。
"""
module PNJLCore

using Base.MathConstants: π
using StaticArrays

const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

const _CONSTANTS_PATH = normpath(joinpath(@__DIR__, "..", "..", "constants", "Constants_PNJL.jl"))
IncludeOnce.include_once!(Main, :Constants_PNJL, _CONSTANTS_PATH)

using Main.Constants_PNJL: pnjl_constants

export PNJLParams
export pnjl_params, as_namedtuple
export calculate_mass_vec, chiral_potential, polyakov_potential
export vacuum_integral_with_cutoff
export cached_nodes, calculate_log_sum
export DEFAULT_MOMENTUM_COUNT, DEFAULT_THETA_COUNT

const _PNJL_INTEGRALS_PATH = normpath(joinpath(@__DIR__, "PNJLIntegrals.jl"))
if !isdefined(@__MODULE__, :PNJLIntegrals)
    include(_PNJL_INTEGRALS_PATH)
end

const _POLYAKOV_EPS = 1e-16
const DEFAULT_MOMENTUM_COUNT = PNJLIntegrals.DEFAULT_MOMENTUM_COUNT
const DEFAULT_THETA_COUNT = PNJLIntegrals.DEFAULT_THETA_COUNT

Base.@kwdef struct PNJLParams
    profile::String = "default"
    physics_profile::String = "default"
    path::Union{Nothing, String} = nothing

    hbarc_MeV_fm::Float64 = 197.327
    alpha_em::Float64 = 1.0 / 137.035999084

    N_color::Int = 3
    N_flavor::Int = 3
    rho0_fm3::Float64 = 0.16

    Λ_inv_fm::Float64
    m_ud0_inv_fm::Float64
    m_s0_inv_fm::Float64
    G_fm2::Float64
    K_fm5::Float64

    T0_inv_fm::Float64
    a0::Float64
    a1::Float64
    a2::Float64
    b3::Float64
    b4::Float64
end

function pnjl_params(constants::NamedTuple)
    return PNJLParams(
        profile=String(get(constants, :profile, "default")),
        physics_profile=String(get(constants, :physics_profile, "default")),
        path=get(constants, :path, nothing),
        hbarc_MeV_fm=Float64(constants.hbarc_MeV_fm),
        alpha_em=Float64(constants.alpha_em),
        N_color=Int(constants.N_color),
        N_flavor=Int(constants.N_flavor),
        rho0_fm3=Float64(constants.rho0_fm3),
        Λ_inv_fm=Float64(constants.Λ_inv_fm),
        m_ud0_inv_fm=Float64(constants.m_ud0_inv_fm),
        m_s0_inv_fm=Float64(constants.m_s0_inv_fm),
        G_fm2=Float64(constants.G_fm2),
        K_fm5=Float64(constants.K_fm5),
        T0_inv_fm=Float64(constants.T0_inv_fm),
        a0=Float64(constants.a0),
        a1=Float64(constants.a1),
        a2=Float64(constants.a2),
        b3=Float64(constants.b3),
        b4=Float64(constants.b4),
    )
end

function pnjl_params(; profile::String=get(ENV, "PNJL_PARAM_PROFILE", "default"), physics_profile::String=get(ENV, "PHYSICS_PARAM_PROFILE", "default"), log_config::Bool=false)
    return pnjl_params(pnjl_constants(profile=profile, physics_profile=physics_profile, log_config=log_config))
end

@inline function as_namedtuple(p::PNJLParams)
    return (
        profile=p.profile,
        physics_profile=p.physics_profile,
        path=p.path,
        hbarc_MeV_fm=p.hbarc_MeV_fm,
        alpha_em=p.alpha_em,
        N_color=p.N_color,
        N_flavor=p.N_flavor,
        rho0_fm3=p.rho0_fm3,
        Λ_inv_fm=p.Λ_inv_fm,
        m_ud0_inv_fm=p.m_ud0_inv_fm,
        m_s0_inv_fm=p.m_s0_inv_fm,
        G_fm2=p.G_fm2,
        K_fm5=p.K_fm5,
        T0_inv_fm=p.T0_inv_fm,
        a0=p.a0,
        a1=p.a1,
        a2=p.a2,
        b3=p.b3,
        b4=p.b4,
    )
end

function calculate_mass_vec(p::PNJLParams, φ::SVector{3, T}) where {T}
    φ_u, φ_d, φ_s = φ
    m_ud0 = convert(T, p.m_ud0_inv_fm)
    m_s0 = convert(T, p.m_s0_inv_fm)
    G = convert(T, p.G_fm2)
    K = convert(T, p.K_fm5)

    return SVector{3, T}(
        m_ud0 - 4 * G * φ_u + 2 * K * φ_d * φ_s,
        m_ud0 - 4 * G * φ_d + 2 * K * φ_u * φ_s,
        m_s0 - 4 * G * φ_s + 2 * K * φ_u * φ_d,
    )
end

@inline function chiral_potential(p::PNJLParams, φ::SVector{3, T}) where {T}
    G = convert(T, p.G_fm2)
    K = convert(T, p.K_fm5)
    return 2 * G * sum(φ .^ 2) - 4 * K * prod(φ)
end

@inline function _safe_log(x; min_val=_POLYAKOV_EPS)
    min_x = one(x) * min_val
    x <= zero(x) && return log(min_x)
    return x < min_x ? log(min_x) : log(x)
end

function polyakov_potential(p::PNJLParams, Φ, Φbar, T_fm)
    TT = promote_type(typeof(Φ), typeof(Φbar), typeof(T_fm))
    ΦT = convert(TT, Φ)
    ΦbarT = convert(TT, Φbar)
    TT_fm = convert(TT, T_fm)

    T0 = convert(TT, p.T0_inv_fm)
    a0 = convert(TT, p.a0)
    a1 = convert(TT, p.a1)
    a2 = convert(TT, p.a2)
    b3 = convert(TT, p.b3)

    T_ratio = T0 / TT_fm
    Ta = a0 + a1 * T_ratio + a2 * T_ratio^2
    Tb = b3 * T_ratio^3
    value = 1 - 6 * ΦbarT * ΦT + 4 * (ΦbarT^3 + ΦT^3) - 3 * (ΦbarT * ΦT)^2
    return TT_fm^4 * (-0.5 * Ta * ΦbarT * ΦT + Tb * _safe_log(value))
end

@inline function vacuum_integral_with_cutoff(mass::T, Λ::T) where {T}
    mass_abs = abs(mass)
    epsilon = one(mass_abs) * 1e-12
    mass_safe = mass_abs + epsilon
    sqrt_term = sqrt(Λ^2 + mass_safe^2)
    poly_part = Λ * sqrt_term * (2 * Λ^2 + mass_safe^2)
    log_term = mass_safe^4 * log((Λ + sqrt_term) / mass_safe)
    return (poly_part - log_term) / (16 * π^2)
end

@inline cached_nodes(p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT) =
    PNJLIntegrals.cached_nodes(p_num, t_num)

@inline calculate_log_sum(args...) = PNJLIntegrals.calculate_log_sum(args...)

end # module PNJLCore
