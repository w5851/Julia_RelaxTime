module TransportCoefficients

"""
在弛豫时间近似（RTA）下计算夸克物质输运系数：
- 剪切粘滞系数 η
- 体粘滞系数 ζ
- 电导率 σ

单位约定（自然单位 ħ=c=k_B=1）：
- 温度/化学势/质量/动量均使用 fm⁻¹
- τ 使用 fm

实现约定：
- τ 与动量 p 无关（直接使用 `RelaxationTime.relaxation_times` 的输出）
- ξ=0：各向同性，仅对动量 p 积分
  - 相空间测度：4π/(2π)³ ∫ p² dp = 1/(2π²) ∫ p² dp
- ξ≠0：分布函数采用 RS 形式，需要完整的角度积分
  - 相空间测度：2π/(2π)³ ∫ p² dp ∫ d(cosθ) = 1/(4π²) ∫ p² dp ∫ d(cosθ)

积分核（含相空间测度 p²）：
- η: p⁶/E² (来自 p² × p⁴/E²)
- σ: p⁴q²/E² (来自 p² × p²q²/E²)
"""

# Include-once helper
const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

# Prefer reuse Main.Constants_PNJL to avoid duplicating the constants module
const _CONSTANTS_PNJL_PATH = normpath(joinpath(@__DIR__, "..", "constants", "Constants_PNJL.jl"))
IncludeOnce.include_once!(Main, :Constants_PNJL, _CONSTANTS_PNJL_PATH)
include("../integration/GaussLegendre.jl")
include("../integration/PhaseSpaceSampling.jl")

# Prefer reuse Main-level distribution modules to avoid duplication
const _QUARK_DISTRIBUTION_PATH = normpath(joinpath(@__DIR__, "..", "models", "pnjl_physics", "QuarkDistribution.jl"))
IncludeOnce.include_once!(Main, :PNJLQuarkDistributions, _QUARK_DISTRIBUTION_PATH)
const _QUARK_DISTRIBUTION_ANISO_PATH = normpath(joinpath(@__DIR__, "QuarkDistribution_Aniso.jl"))
IncludeOnce.include_once!(Main, :PNJLQuarkDistributions_Aniso, _QUARK_DISTRIBUTION_ANISO_PATH)

using Base.MathConstants: π

using Main.Constants_PNJL: N_color, α
using .GaussLegendre: gauleg, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS, DEFAULT_COSΘ_NODES, DEFAULT_COSΘ_WEIGHTS
using .PhaseSpaceSampling: p_nodes_weights, cos_nodes_weights, integrate_p, integrate_p_cos
using Main.PNJLQuarkDistributions: quark_distribution, antiquark_distribution
using Main.PNJLQuarkDistributions_Aniso: quark_distribution_aniso, antiquark_distribution_aniso

# Ensure shared parameter types are loaded at Main-level
const _PARAMETER_TYPES_PATH = normpath(joinpath(@__DIR__, "..", "types", "ParameterTypes.jl"))
IncludeOnce.include_once!(Main, :ParameterTypes, _PARAMETER_TYPES_PATH)

using Main.ParameterTypes: QuarkParams, ThermoParams

export shear_viscosity, bulk_viscosity, bulk_viscosity_isentropic, electric_conductivity, transport_coefficients

export TransportIntegrationConfig, QuarkParams, ThermoParams, TransportPhysicsConfig, TransportRequest
export default_transport_provider

const TWO_PI = 2.0 * π
const ENERGY_FLOOR = sqrt(eps(Float64))

const _SPECIES_ALL = (:u, :d, :s, :ubar, :dbar, :sbar)
const _FLAVORS = (:u, :d, :s)

# 自然单位制中的元电荷: e = sqrt(4πα)
const e_natural = sqrt(4.0 * π * α)

const _TRANSPORT_INTEGRATION_KEYS = (
    :p_nodes,
    :p_max,
    :p_grid,
    :p_w,
    :cos_nodes,
    :cos_grid,
    :cos_w,
)

struct TransportIntegrationConfig
    p_nodes::Int
    p_max::Float64
    p_grid::Union{Nothing,Vector{Float64}}
    p_w::Union{Nothing,Vector{Float64}}
    cos_nodes::Int
    cos_grid::Union{Nothing,Vector{Float64}}
    cos_w::Union{Nothing,Vector{Float64}}
end

function TransportIntegrationConfig(;
    p_nodes::Int=length(DEFAULT_MOMENTUM_NODES),
    p_max::Float64=10.0,
    p_grid::Union{Nothing,Vector{Float64}}=nothing,
    p_w::Union{Nothing,Vector{Float64}}=nothing,
    cos_nodes::Int=length(DEFAULT_COSΘ_NODES),
    cos_grid::Union{Nothing,Vector{Float64}}=nothing,
    cos_w::Union{Nothing,Vector{Float64}}=nothing,
)
    p_nodes > 0 || error("TransportIntegrationConfig: p_nodes must be > 0")
    cos_nodes > 0 || error("TransportIntegrationConfig: cos_nodes must be > 0")
    isfinite(p_max) && p_max > 0.0 || error("TransportIntegrationConfig: p_max must be finite and > 0")

    if (p_grid === nothing) != (p_w === nothing)
        error("TransportIntegrationConfig: p_grid and p_w must be provided together")
    end
    if p_grid !== nothing && length(p_grid) != length(p_w)
        error("TransportIntegrationConfig: p_grid and p_w length mismatch")
    end
    if p_grid !== nothing
        isempty(p_grid) && error("TransportIntegrationConfig: p_grid must not be empty")
        all(isfinite, p_grid) || error("TransportIntegrationConfig: p_grid must be finite")
        all(isfinite, p_w) || error("TransportIntegrationConfig: p_w must be finite")
    end
    if (cos_grid === nothing) != (cos_w === nothing)
        error("TransportIntegrationConfig: cos_grid and cos_w must be provided together")
    end
    if cos_grid !== nothing && length(cos_grid) != length(cos_w)
        error("TransportIntegrationConfig: cos_grid and cos_w length mismatch")
    end
    if cos_grid !== nothing
        isempty(cos_grid) && error("TransportIntegrationConfig: cos_grid must not be empty")
        all(isfinite, cos_grid) || error("TransportIntegrationConfig: cos_grid must be finite")
        all(isfinite, cos_w) || error("TransportIntegrationConfig: cos_w must be finite")
    end
    return TransportIntegrationConfig(p_nodes, p_max, p_grid, p_w, cos_nodes, cos_grid, cos_w)
end

const DEFAULT_TRANSPORT_CONFIG = TransportIntegrationConfig()

@inline function _merge_transport_integration_config(
    base::TransportIntegrationConfig,
    overrides::Base.Pairs,
)::TransportIntegrationConfig
    unknown = Symbol[]
    for k in keys(overrides)
        k in _TRANSPORT_INTEGRATION_KEYS || push!(unknown, k)
    end
    isempty(unknown) || error(
        "Unknown integration keyword(s): $(unknown). Allowed keys: $(_TRANSPORT_INTEGRATION_KEYS)"
    )

    return TransportIntegrationConfig(
        p_nodes=get(overrides, :p_nodes, base.p_nodes),
        p_max=get(overrides, :p_max, base.p_max),
        p_grid=get(overrides, :p_grid, base.p_grid),
        p_w=get(overrides, :p_w, base.p_w),
        cos_nodes=get(overrides, :cos_nodes, base.cos_nodes),
        cos_grid=get(overrides, :cos_grid, base.cos_grid),
        cos_w=get(overrides, :cos_w, base.cos_w),
    )
end

@inline energy_from_p(p::Float64, m::Float64) = sqrt(p * p + m * m)

@inline function energy_from_p_aniso(p::Float64, m::Float64, ξ::Float64, cosθ::Float64)
    return sqrt(p * p + m * m + ξ * (p * cosθ)^2)
end

const DEFAULT_TRANSPORT_PROVIDER = (
    energy_from_p=energy_from_p,
    energy_from_p_aniso=energy_from_p_aniso,
    quark_distribution=quark_distribution,
    antiquark_distribution=antiquark_distribution,
    quark_distribution_aniso=quark_distribution_aniso,
    antiquark_distribution_aniso=antiquark_distribution_aniso,
    prefer_energy_aniso=true,
)

@inline default_transport_provider() = DEFAULT_TRANSPORT_PROVIDER

@inline function degeneracy_default()::Float64
    return 2.0 * Float64(N_color)
end

"""
    default_charges()

返回夸克电荷（自然单位制）。

在自然单位制中，元电荷 e = √(4πα)，其中 α ≈ 1/137 是精细结构常数。
夸克电荷为：
- u夸克: q_u = (2/3)e
- d夸克: q_d = (-1/3)e  
- s夸克: q_s = (-1/3)e
"""
@inline function default_charges()
    return (u = (2.0 / 3.0) * e_natural, d = (-1.0 / 3.0) * e_natural, s = (-1.0 / 3.0) * e_natural)
end

struct TransportPhysicsConfig
    charges::NamedTuple
    degeneracy::Float64
end

function TransportPhysicsConfig(;
    charges::NamedTuple=default_charges(),
    degeneracy::Float64=degeneracy_default(),
)
    return TransportPhysicsConfig(charges, degeneracy)
end

struct TransportRequest
    quark::QuarkParams
    thermo::ThermoParams
    tau::NamedTuple
    physics::TransportPhysicsConfig
    integration::TransportIntegrationConfig
end

function TransportRequest(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    tau::NamedTuple,
    charges::NamedTuple=default_charges(),
    degeneracy::Float64=degeneracy_default(),
    integration::TransportIntegrationConfig=DEFAULT_TRANSPORT_CONFIG,
)
    return TransportRequest(
        QuarkParams(quark_params),
        ThermoParams(thermo_params),
        tau,
        TransportPhysicsConfig(charges=charges, degeneracy=degeneracy),
        integration,
    )
end

function TransportRequest(
    quark::QuarkParams,
    thermo::ThermoParams;
    tau::NamedTuple,
    physics::TransportPhysicsConfig=TransportPhysicsConfig(),
    integration::TransportIntegrationConfig=DEFAULT_TRANSPORT_CONFIG,
)
    return TransportRequest(quark, thermo, tau, physics, integration)
end

@inline _qp_view(q::QuarkParams) = (m=q.m, μ=q.μ)
@inline _tp_view(t::ThermoParams) = (T=t.T, Φ=t.Φ, Φbar=t.Φbar, ξ=t.ξ)

function shear_viscosity(
    quark::QuarkParams,
    thermo::ThermoParams;
    kwargs...
)::Float64
    return shear_viscosity(_qp_view(quark), _tp_view(thermo); kwargs...)
end

function shear_viscosity(
    quark::QuarkParams,
    thermo::ThermoParams,
    config::TransportIntegrationConfig;
    kwargs...
)::Float64
    return shear_viscosity(_qp_view(quark), _tp_view(thermo), config; kwargs...)
end

function electric_conductivity(
    quark::QuarkParams,
    thermo::ThermoParams;
    kwargs...
)::Float64
    return electric_conductivity(_qp_view(quark), _tp_view(thermo); kwargs...)
end

function electric_conductivity(
    quark::QuarkParams,
    thermo::ThermoParams,
    config::TransportIntegrationConfig;
    kwargs...
)::Float64
    return electric_conductivity(_qp_view(quark), _tp_view(thermo), config; kwargs...)
end

function bulk_viscosity_isentropic(
    quark::QuarkParams,
    thermo::ThermoParams;
    kwargs...
)::Float64
    return bulk_viscosity_isentropic(_qp_view(quark), _tp_view(thermo); kwargs...)
end

function bulk_viscosity_isentropic(
    quark::QuarkParams,
    thermo::ThermoParams,
    config::TransportIntegrationConfig;
    kwargs...
)::Float64
    return bulk_viscosity_isentropic(_qp_view(quark), _tp_view(thermo), config; kwargs...)
end

"""Unified bulk viscosity entry. Currently supports `formula=:isentropic`."""
function bulk_viscosity(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    formula::Symbol=:isentropic,
    bulk_coeffs_isentropic::Union{Nothing,NamedTuple}=nothing,
    bulk_coeffs::Union{Nothing,NamedTuple}=nothing,
    kwargs...
)::Float64
    formula === :isentropic || error("Unsupported bulk viscosity formula: $formula. Supported formulas: (:isentropic,)")
    coeffs = something(bulk_coeffs_isentropic, bulk_coeffs)
    coeffs === nothing && error("bulk_viscosity requires bulk_coeffs_isentropic (or legacy alias bulk_coeffs) when formula=:isentropic")
    return bulk_viscosity_isentropic(quark_params, thermo_params; bulk_coeffs_isentropic=coeffs, kwargs...)
end

function bulk_viscosity(
    quark::QuarkParams,
    thermo::ThermoParams;
    kwargs...
)::Float64
    return bulk_viscosity(_qp_view(quark), _tp_view(thermo); kwargs...)
end

function bulk_viscosity(
    req::TransportRequest;
    formula::Symbol=:isentropic,
    kwargs...
)::Float64
    quark_params = (m=req.quark.m, μ=req.quark.μ)
    thermo_params = (T=req.thermo.T, Φ=req.thermo.Φ, Φbar=req.thermo.Φbar, ξ=req.thermo.ξ)
    return bulk_viscosity(quark_params, thermo_params; formula=formula, tau=req.tau, kwargs...)
end

function transport_coefficients(
    quark::QuarkParams,
    thermo::ThermoParams;
    kwargs...
)::NamedTuple
    return transport_coefficients(_qp_view(quark), _tp_view(thermo); kwargs...)
end

function transport_coefficients(
    quark::QuarkParams,
    thermo::ThermoParams,
    config::TransportIntegrationConfig;
    kwargs...
)::NamedTuple
    return transport_coefficients(_qp_view(quark), _tp_view(thermo), config; kwargs...)
end

"""Recommended entry: pass a `TransportRequest` and optionally override integration fields by keywords."""
function shear_viscosity(
    req::TransportRequest;
    tau::NamedTuple=req.tau,
    degeneracy::Float64=req.physics.degeneracy,
    integration::TransportIntegrationConfig=req.integration,
    provider=DEFAULT_TRANSPORT_PROVIDER,
    kwargs...
)::Float64
    quark_params = (m=req.quark.m, μ=req.quark.μ)
    thermo_params = (T=req.thermo.T, Φ=req.thermo.Φ, Φbar=req.thermo.Φbar, ξ=req.thermo.ξ)
    return shear_viscosity(
        quark_params,
        thermo_params;
        tau=tau,
        degeneracy=degeneracy,
        config=integration,
        provider=provider,
        kwargs...
    )
end

"""Recommended entry: pass a `TransportRequest` and optionally override integration fields by keywords."""
function electric_conductivity(
    req::TransportRequest;
    tau::NamedTuple=req.tau,
    charges::NamedTuple=req.physics.charges,
    degeneracy::Float64=req.physics.degeneracy,
    integration::TransportIntegrationConfig=req.integration,
    provider=DEFAULT_TRANSPORT_PROVIDER,
    kwargs...
)::Float64
    quark_params = (m=req.quark.m, μ=req.quark.μ)
    thermo_params = (T=req.thermo.T, Φ=req.thermo.Φ, Φbar=req.thermo.Φbar, ξ=req.thermo.ξ)
    return electric_conductivity(
        quark_params,
        thermo_params;
        tau=tau,
        charges=charges,
        degeneracy=degeneracy,
        config=integration,
        provider=provider,
        kwargs...
    )
end

"""Recommended entry: pass a `TransportRequest` and optionally override integration fields by keywords."""
function bulk_viscosity_isentropic(
    req::TransportRequest;
    tau::NamedTuple=req.tau,
    bulk_coeffs_isentropic::NamedTuple,
    integration::TransportIntegrationConfig=req.integration,
    provider=DEFAULT_TRANSPORT_PROVIDER,
    kwargs...
)::Float64
    quark_params = (m=req.quark.m, μ=req.quark.μ)
    thermo_params = (T=req.thermo.T, Φ=req.thermo.Φ, Φbar=req.thermo.Φbar, ξ=req.thermo.ξ)
    return bulk_viscosity_isentropic(
        quark_params,
        thermo_params;
        tau=tau,
        bulk_coeffs_isentropic=bulk_coeffs_isentropic,
        config=integration,
        provider=provider,
        kwargs...
    )
end

"""Recommended entry: pass a `TransportRequest` and optionally override integration fields by keywords."""
function transport_coefficients(
    req::TransportRequest;
    tau::NamedTuple=req.tau,
    bulk_coeffs::Union{Nothing,NamedTuple}=nothing,
    charges::NamedTuple=req.physics.charges,
    degeneracy::Float64=req.physics.degeneracy,
    integration::TransportIntegrationConfig=req.integration,
    provider=DEFAULT_TRANSPORT_PROVIDER,
    kwargs...
)::NamedTuple
    quark_params = (m=req.quark.m, μ=req.quark.μ)
    thermo_params = (T=req.thermo.T, Φ=req.thermo.Φ, Φbar=req.thermo.Φbar, ξ=req.thermo.ξ)
    return transport_coefficients(
        quark_params,
        thermo_params;
        tau=tau,
        bulk_coeffs=bulk_coeffs,
        charges=charges,
        degeneracy=degeneracy,
        config=integration,
        provider=provider,
        kwargs...
    )
end

@inline function fermi_factor(f::Float64)
    return clamp(f * (1.0 - f), 0.0, 0.25)
end

@inline function _is_quark_species(species::Symbol)::Bool
    return species === :u || species === :d || species === :s
end

@inline function _flavor_index(species::Symbol)::Int
    if species === :u || species === :ubar
        return 1
    elseif species === :d || species === :dbar
        return 2
    elseif species === :s || species === :sbar
        return 3
    end
    error("Unknown species: $species")
end

@inline function _anti_species(flavor::Symbol)::Symbol
    if flavor === :u
        return :ubar
    elseif flavor === :d
        return :dbar
    elseif flavor === :s
        return :sbar
    end
    error("Unknown flavor: $flavor")
end

@inline function _q2_for_species(species::Symbol, charges::NamedTuple)::Float64
    idx = _flavor_index(species)
    if idx == 1
        return charges.u^2
    elseif idx == 2
        return charges.d^2
    end
    return charges.s^2
end

@inline function _provider_caps(provider)
    return (
        has_energy_from_p_aniso=hasproperty(provider, :energy_from_p_aniso),
        has_quark_distribution_aniso=hasproperty(provider, :quark_distribution_aniso),
        has_antiquark_distribution_aniso=hasproperty(provider, :antiquark_distribution_aniso),
        has_mass_for_species=hasproperty(provider, :mass_for_species),
        has_mu_for_species=hasproperty(provider, :mu_for_species),
        prefer_energy_aniso=(hasproperty(provider, :prefer_energy_aniso) && provider.prefer_energy_aniso === true),
    )
end

@inline function _energy_for_kernel(provider, caps::NamedTuple, p::Float64, m::Float64, ξ::Float64, cosθ::Float64)::Float64
    raw_E = if ξ == 0.0 || !caps.has_energy_from_p_aniso
        provider.energy_from_p(p, m)
    else
        provider.energy_from_p_aniso(p, m, ξ, cosθ)
    end
    return max(raw_E, ENERGY_FLOOR)
end

@inline function distribution_for_species_from_E(
    provider,
    caps::NamedTuple,
    species::Symbol,
    E::Float64,
    p::Float64,
    m::Float64,
    μ::Float64,
    T::Float64,
    Φ::Float64,
    Φbar::Float64,
    ξ::Float64,
    cosθ::Float64
)
    if ξ == 0.0
        if _is_quark_species(species)
            return provider.quark_distribution(E, μ, T, Φ, Φbar)
        else
            return provider.antiquark_distribution(E, μ, T, Φ, Φbar)
        end
    else
        is_quark = _is_quark_species(species)

        has_energy_passthrough = caps.has_energy_from_p_aniso
        has_aniso_distribution = is_quark ? caps.has_quark_distribution_aniso : caps.has_antiquark_distribution_aniso

        if has_energy_passthrough && (caps.prefer_energy_aniso || !has_aniso_distribution)
            return is_quark ? provider.quark_distribution(E, μ, T, Φ, Φbar) : provider.antiquark_distribution(E, μ, T, Φ, Φbar)
        end

        if has_aniso_distribution
            return is_quark ? provider.quark_distribution_aniso(p, m, μ, T, Φ, Φbar, ξ, cosθ) : provider.antiquark_distribution_aniso(p, m, μ, T, Φ, Φbar, ξ, cosθ)
        end

        error(
            "provider does not support anisotropic distribution: missing $(is_quark ? :quark_distribution_aniso : :antiquark_distribution_aniso) and no energy passthrough available (need :energy_from_p_aniso + :quark_distribution/:antiquark_distribution)"
        )
    end
end

@inline function distribution_for_species_from_E(
    provider,
    species::Symbol,
    E::Float64,
    p::Float64,
    m::Float64,
    μ::Float64,
    T::Float64,
    Φ::Float64,
    Φbar::Float64,
    ξ::Float64,
    cosθ::Float64
)
    caps = _provider_caps(provider)
    return distribution_for_species_from_E(provider, caps, species, E, p, m, μ, T, Φ, Φbar, ξ, cosθ)
end

@inline function distribution_for_species(
    provider,
    species::Symbol,
    p::Float64,
    m::Float64,
    μ::Float64,
    T::Float64,
    Φ::Float64,
    Φbar::Float64,
    ξ::Float64,
    cosθ::Float64
)
    caps = _provider_caps(provider)
    E = _energy_for_kernel(provider, caps, p, m, ξ, cosθ)
    return distribution_for_species_from_E(provider, species, E, p, m, μ, T, Φ, Φbar, ξ, cosθ)
end

@inline function _energy_for_kernel(provider, p::Float64, m::Float64, ξ::Float64, cosθ::Float64)::Float64
    caps = _provider_caps(provider)
    return _energy_for_kernel(provider, caps, p, m, ξ, cosθ)
end

@inline function mass_for_species(species::Symbol, quark_params::Union{NamedTuple,QuarkParams})::Float64
    idx = _flavor_index(species)
    if idx == 1
        return quark_params.m.u
    elseif idx == 2
        return quark_params.m.d
    else
        return quark_params.m.s
    end
end

@inline function mu_for_species(species::Symbol, quark_params::Union{NamedTuple,QuarkParams})::Float64
    idx = _flavor_index(species)
    if idx == 1
        return quark_params.μ.u
    elseif idx == 2
        return quark_params.μ.d
    elseif idx == 3
        return quark_params.μ.s
    end
    error("Unknown species: $species")
end

@inline function _mass_for_species(provider, species::Symbol, quark_params, thermo_params)::Float64
    if hasproperty(provider, :mass_for_species)
        return provider.mass_for_species(species, quark_params, thermo_params)
    end
    return mass_for_species(species, quark_params)
end

@inline function _mass_for_species(provider, caps::NamedTuple, species::Symbol, quark_params, thermo_params)::Float64
    if caps.has_mass_for_species
        return provider.mass_for_species(species, quark_params, thermo_params)
    end
    return mass_for_species(species, quark_params)
end

@inline function _mu_for_species(provider, species::Symbol, quark_params, thermo_params)::Float64
    if hasproperty(provider, :mu_for_species)
        return provider.mu_for_species(species, quark_params, thermo_params)
    end
    return mu_for_species(species, quark_params)
end

@inline function _mu_for_species(provider, caps::NamedTuple, species::Symbol, quark_params, thermo_params)::Float64
    if caps.has_mu_for_species
        return provider.mu_for_species(species, quark_params, thermo_params)
    end
    return mu_for_species(species, quark_params)
end

@inline function tau_for_species(species::Symbol, tau::NamedTuple)::Float64
    hasproperty(tau, species) || error("tau is missing :$species")
    return getproperty(tau, species)
end

@inline function _species_transport_state(
    provider,
    caps::NamedTuple,
    sp::Symbol,
    p::Float64,
    c::Float64,
    ξ::Float64,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    tau::NamedTuple,
)::NTuple{3,Float64}
    T = thermo_params.T
    Φ = thermo_params.Φ
    Φbar = thermo_params.Φbar

    m = _mass_for_species(provider, caps, sp, quark_params, thermo_params)
    μ = _mu_for_species(provider, caps, sp, quark_params, thermo_params)
    E = _energy_for_kernel(provider, caps, p, m, ξ, c)
    f = distribution_for_species_from_E(provider, caps, sp, E, p, m, μ, T, Φ, Φbar, ξ, c)
    ff = fermi_factor(f)
    τ = tau_for_species(sp, tau)
    return (E, ff, τ)
end

@inline _p_nodes_weights(p_nodes::Int, p_max::Float64, p_grid, p_w) = p_nodes_weights(p_nodes, p_max, p_grid, p_w)
@inline _cos_nodes_weights(cos_nodes::Int, cos_grid, cos_w) = cos_nodes_weights(cos_nodes, cos_grid, cos_w)

@inline function _effective_transport_config(
    config::Union{Nothing,TransportIntegrationConfig},
    kwargs::Base.Pairs,
)::TransportIntegrationConfig
    base_config = something(config, DEFAULT_TRANSPORT_CONFIG)
    return _merge_transport_integration_config(base_config, kwargs)
end

@inline function _validate_transport_inputs(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    tau::NamedTuple,
    config::TransportIntegrationConfig;
    provider=nothing,
    charges::Union{Nothing,NamedTuple}=nothing,
    bulk_coeffs_isentropic::Union{Nothing,NamedTuple}=nothing,
)
    T = thermo_params.T
    isfinite(T) && T > 0.0 || error("thermo_params.T must be finite and > 0")
    isfinite(thermo_params.Φ) || error("thermo_params.Φ must be finite")
    isfinite(thermo_params.Φbar) || error("thermo_params.Φbar must be finite")

    for sp in _SPECIES_ALL
        τ = tau_for_species(sp, tau)
        isfinite(τ) || error("tau.:$sp must be finite")
        τ >= 0.0 || error("tau.:$sp must be >= 0")
    end

    check_mass = provider === nothing || !hasproperty(provider, :mass_for_species)
    check_mu = provider === nothing || !hasproperty(provider, :mu_for_species)

    for fl in _FLAVORS
        m = getproperty(quark_params.m, fl)
        μ = getproperty(quark_params.μ, fl)
        if check_mass
            isfinite(m) || error("quark_params.m.$fl must be finite")
            m >= 0.0 || error("quark_params.m.$fl must be >= 0")
        end
        if check_mu
            isfinite(μ) || error("quark_params.μ.$fl must be finite")
        end
    end

    config.p_nodes > 0 || error("integration p_nodes must be > 0")
    config.cos_nodes > 0 || error("integration cos_nodes must be > 0")
    isfinite(config.p_max) && config.p_max > 0.0 || error("integration p_max must be finite and > 0")

    if charges !== nothing
        isfinite(charges.u) || error("charges.u must be finite")
        isfinite(charges.d) || error("charges.d must be finite")
        isfinite(charges.s) || error("charges.s must be finite")
    end

    if bulk_coeffs_isentropic !== nothing
        isfinite(bulk_coeffs_isentropic.v_n_sq) || error("bulk_coeffs_isentropic.v_n_sq must be finite")
        isfinite(bulk_coeffs_isentropic.dμB_dT_sigma) || error("bulk_coeffs_isentropic.dμB_dT_sigma must be finite")

        length(bulk_coeffs_isentropic.masses) == 3 || error("bulk_coeffs_isentropic.masses must have length 3")
        length(bulk_coeffs_isentropic.dM_dT) == 3 || error("bulk_coeffs_isentropic.dM_dT must have length 3")
        length(bulk_coeffs_isentropic.dM_dμB) == 3 || error("bulk_coeffs_isentropic.dM_dμB must have length 3")

        all(isfinite, bulk_coeffs_isentropic.masses) || error("bulk_coeffs_isentropic.masses must be finite")
        all(m -> m >= 0.0, bulk_coeffs_isentropic.masses) || error("bulk_coeffs_isentropic.masses must be >= 0")
        all(isfinite, bulk_coeffs_isentropic.dM_dT) || error("bulk_coeffs_isentropic.dM_dT must be finite")
        all(isfinite, bulk_coeffs_isentropic.dM_dμB) || error("bulk_coeffs_isentropic.dM_dμB must be finite")
    end

    return nothing
end

@inline function _integrate_species_sum(
    term_for_species,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    tau::NamedTuple,
    provider,
    config::TransportIntegrationConfig,
)::Float64
    ξ = get(thermo_params, :ξ, 0.0)
    caps = _provider_caps(provider)

    nodes_p, weights_p = _p_quadrature(config)
    pref_measure_iso = 1.0 / (2.0 * π^2)
    pref_measure_aniso = 1.0 / (4.0 * π^2)

    if ξ == 0.0
        acc = integrate_p(nodes_p, weights_p) do p
            inner = 0.0
            for sp in _SPECIES_ALL
                E, ff, τ = _species_transport_state(provider, caps, sp, p, 0.0, 0.0, quark_params, thermo_params, tau)
                inner += term_for_species(p, sp, E, ff, τ)
            end
            return inner
        end
        return pref_measure_iso * acc
    end

    nodes_cos, weights_cos = _cos_quadrature(config)
    acc = integrate_p_cos(nodes_p, weights_p, nodes_cos, weights_cos) do p, c
        inner = 0.0
        for sp in _SPECIES_ALL
            E, ff, τ = _species_transport_state(provider, caps, sp, p, c, ξ, quark_params, thermo_params, tau)
            inner += term_for_species(p, sp, E, ff, τ)
        end
        return inner
    end
    return pref_measure_aniso * acc
end

@inline function _integrate_species_sum(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    tau::NamedTuple,
    provider,
    config::TransportIntegrationConfig,
    term_for_species,
)::Float64
    return _integrate_species_sum(term_for_species, quark_params, thermo_params, tau, provider, config)
end

@inline function _p_quadrature(effective_config::TransportIntegrationConfig)
    return _p_nodes_weights(
        effective_config.p_nodes,
        effective_config.p_max,
        effective_config.p_grid,
        effective_config.p_w,
    )
end

@inline function _cos_quadrature(effective_config::TransportIntegrationConfig)
    return _cos_nodes_weights(
        effective_config.cos_nodes,
        effective_config.cos_grid,
        effective_config.cos_w,
    )
end

"""
    shear_viscosity(quark_params, thermo_params; tau, ...)

剪切粘滞系数 η（RTA）。

积分公式：
- 各向同性 (ξ=0): η = (1/15T) × (1/2π²) × ∫ dp p⁶/E² × g × τ × f(1-f)
- 各向异性 (ξ≠0): η = (1/15T) × (1/4π²) × ∫ dp ∫ d(cosθ) p⁶/E² × g × τ × f(1-f)
"""
function shear_viscosity(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    tau::NamedTuple,
    degeneracy::Float64=degeneracy_default(),
    provider=DEFAULT_TRANSPORT_PROVIDER,
    config::Union{Nothing,TransportIntegrationConfig}=nothing,
    kwargs...
)::Float64
    T = thermo_params.T

    effective_config = _effective_transport_config(config, kwargs)
    _validate_transport_inputs(quark_params, thermo_params, tau, effective_config; provider=provider)

    integral = _integrate_species_sum(quark_params, thermo_params, tau, provider, effective_config) do p, sp, E, ff, τ
        p2 = p * p
        p6 = p2 * p2 * p2
        return p6 / (E * E) * (degeneracy * τ * ff)
    end

    return (1.0 / (15.0 * T)) * integral
end

"""Convenience overload: pass `config` as a positional argument."""
function shear_viscosity(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    config::TransportIntegrationConfig;
    tau::NamedTuple,
    degeneracy::Float64=degeneracy_default(),
    provider=DEFAULT_TRANSPORT_PROVIDER,
    kwargs...
)::Float64
    return shear_viscosity(
        quark_params,
        thermo_params;
        tau=tau,
        config=config,
        degeneracy=degeneracy,
        provider=provider,
        kwargs...
    )
end

"""
    electric_conductivity(quark_params, thermo_params; tau, charges, ...)

电导率 σ（RTA）。

积分公式：
- 各向同性 (ξ=0): σ = (1/3T) × (1/2π²) × ∫ dp p⁴q²/E² × g × τ × f(1-f)
- 各向异性 (ξ≠0): σ = (1/3T) × (1/4π²) × ∫ dp ∫ d(cosθ) p⁴q²/E² × g × τ × f(1-f)
"""
function electric_conductivity(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    tau::NamedTuple,
    charges::NamedTuple=default_charges(),
    degeneracy::Float64=degeneracy_default(),
    provider=DEFAULT_TRANSPORT_PROVIDER,
    config::Union{Nothing,TransportIntegrationConfig}=nothing,
    kwargs...
)::Float64
    T = thermo_params.T

    effective_config = _effective_transport_config(config, kwargs)
    _validate_transport_inputs(quark_params, thermo_params, tau, effective_config; provider=provider, charges=charges)

    integral = _integrate_species_sum(quark_params, thermo_params, tau, provider, effective_config) do p, sp, E, ff, τ
        p2 = p * p
        p4 = p2 * p2
        return p4 * _q2_for_species(sp, charges) / (E * E) * (degeneracy * τ * ff)
    end

    return (1.0 / (3.0 * T)) * integral
end

"""Convenience overload: pass `config` as a positional argument."""
function electric_conductivity(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    config::TransportIntegrationConfig;
    tau::NamedTuple,
    charges::NamedTuple=default_charges(),
    degeneracy::Float64=degeneracy_default(),
    provider=DEFAULT_TRANSPORT_PROVIDER,
    kwargs...
)::Float64
    return electric_conductivity(
        quark_params,
        thermo_params;
        tau=tau,
        config=config,
        charges=charges,
        degeneracy=degeneracy,
        provider=provider,
        kwargs...
    )
end

"""
    bulk_viscosity_isentropic(quark_params, thermo_params; tau, bulk_coeffs_isentropic, ...)

体粘滞系数 ζ（RTA，等熵声速形式）。

使用公式 A26（与 Fortran 代码一致）：
    ζ = (N_c)/(9π²T) Σ_a ∫dp (p²/E²)·[τ_a·f_a·(1-f_a)·B_a² + τ_ā·f_ā·(1-f_ā)·B_ā²]

其中 B = p² + 3·v_n²·T²·E·∂[(E∓μ)/T]/∂T|_σ

输入 `bulk_coeffs_isentropic` 应使用 `PNJL.ThermoDerivatives.bulk_viscosity_coefficients` 的返回值。
"""
function bulk_viscosity_isentropic(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    tau::NamedTuple,
    bulk_coeffs_isentropic::NamedTuple,
    provider=DEFAULT_TRANSPORT_PROVIDER,
    config::Union{Nothing,TransportIntegrationConfig}=nothing,
    kwargs...
)::Float64
    T = thermo_params.T
    Φ = thermo_params.Φ
    Φbar = thermo_params.Φbar
    ξ = get(thermo_params, :ξ, 0.0)

    # 从 bulk_coeffs_isentropic 提取系数
    v_n_sq = bulk_coeffs_isentropic.v_n_sq
    dμB_dT_sigma = bulk_coeffs_isentropic.dμB_dT_sigma
    masses = bulk_coeffs_isentropic.masses
    dM_dT = bulk_coeffs_isentropic.dM_dT
    dM_dμB = bulk_coeffs_isentropic.dM_dμB

    effective_config = _effective_transport_config(config, kwargs)
    caps = _provider_caps(provider)
    _validate_transport_inputs(
        quark_params,
        thermo_params,
        tau,
        effective_config;
        provider=provider,
        bulk_coeffs_isentropic=bulk_coeffs_isentropic,
    )
    nodes_p, weights_p = _p_quadrature(effective_config)

    # 系数：N_c / (9π²T)
    # 注意：Fortran 公式中的积分是直接对 p 积分，没有额外的相空间测度因子
    # 所以这里不需要乘 1/(2π²)
    prefactor = Float64(N_color) / (9.0 * π^2 * T)

    flavors = _FLAVORS

    # 计算 B 项
    function compute_B(p::Float64, m::Float64, μ::Float64, dM_dT_val::Float64, dM_dμB_val::Float64, is_antiquark::Bool, ξ::Float64, cosθ::Float64)
        E = _energy_for_kernel(provider, caps, p, m, ξ, cosθ)
        
        # 能量导数（对重子化学势）
        dE_dT = (m / E) * dM_dT_val
        dE_dμB = (m / E) * dM_dμB_val
        
        # 夸克重子数 b_q = 1/3
        b_q = 1.0 / 3.0
        
        if is_antiquark
            # 反夸克：x = E + μ_q
            dx_dT_sigma = dE_dT + (dE_dμB + b_q) * dμB_dT_sigma
            x = E + μ
        else
            # 夸克：x = E - μ_q
            dx_dT_sigma = dE_dT + (dE_dμB - b_q) * dμB_dT_sigma
            x = E - μ
        end
        
        # ∂[x/T]/∂T|_σ
        dxt_dT_sigma = dx_dT_sigma / T - x / T^2
        
        # B = p² + 3·v_n²·T²·E·∂[x/T]/∂T|_σ
        return p * p + 3.0 * v_n_sq * T^2 * E * dxt_dT_sigma
    end

    function one_flavor_pair_contrib(flavor::Symbol, p::Float64, cosθ::Float64)
        sp_q = flavor
        sp_aq = _anti_species(flavor)

        idx = _flavor_index(sp_q)
        m = masses[idx]
        μ = _mu_for_species(provider, caps, sp_q, quark_params, thermo_params)
        E = _energy_for_kernel(provider, caps, p, m, ξ, cosθ)

        f_q = distribution_for_species_from_E(provider, caps, sp_q, E, p, m, μ, T, Φ, Φbar, ξ, cosθ)
        f_aq = distribution_for_species_from_E(provider, caps, sp_aq, E, p, m, μ, T, Φ, Φbar, ξ, cosθ)

        ff_q = fermi_factor(f_q)
        ff_aq = fermi_factor(f_aq)

        τ_q = tau_for_species(sp_q, tau)
        τ_aq = tau_for_species(sp_aq, tau)

        B_q = compute_B(p, m, μ, dM_dT[idx], dM_dμB[idx], false, ξ, cosθ)
        B_aq = compute_B(p, m, μ, dM_dT[idx], dM_dμB[idx], true, ξ, cosθ)

        # 积分核：(p²/E²) × τ × f(1-f) × B²
        kernel_q = (p * p / (E * E)) * τ_q * ff_q * B_q^2
        kernel_aq = (p * p / (E * E)) * τ_aq * ff_aq * B_aq^2

        return kernel_q + kernel_aq
    end

    if ξ == 0.0
        acc = integrate_p(nodes_p, weights_p) do p
            inner = 0.0
            for fl in flavors
                inner += one_flavor_pair_contrib(fl, p, 0.0)
            end
            return inner
        end
        return prefactor * acc
    else
        nodes_cos, weights_cos = _cos_quadrature(effective_config)
        acc = integrate_p_cos(nodes_p, weights_p, nodes_cos, weights_cos) do p, c
            inner = 0.0
            for fl in flavors
                inner += one_flavor_pair_contrib(fl, p, c)
            end
            return inner
        end
        # 各向异性情况下需要除以角度积分的归一化因子 2
        return prefactor * acc / 2.0
    end
end

"""Convenience overload: pass `config` as a positional argument."""
function bulk_viscosity_isentropic(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    config::TransportIntegrationConfig;
    tau::NamedTuple,
    bulk_coeffs_isentropic::NamedTuple,
    provider=DEFAULT_TRANSPORT_PROVIDER,
    kwargs...
)::Float64
    return bulk_viscosity_isentropic(
        quark_params,
        thermo_params;
        tau=tau,
        bulk_coeffs_isentropic=bulk_coeffs_isentropic,
        config=config,
        provider=provider,
        kwargs...
    )
end

"""
    transport_coefficients(quark_params, thermo_params; tau, bulk_coeffs, charges, ...)

一次性计算 (η, ζ, σ)。
"""
function transport_coefficients(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    tau::NamedTuple,
    bulk_coeffs::Union{Nothing,NamedTuple}=nothing,
    config::Union{Nothing,TransportIntegrationConfig}=nothing,
    charges::NamedTuple=default_charges(),
    degeneracy::Float64=degeneracy_default(),
    provider=DEFAULT_TRANSPORT_PROVIDER,
    kwargs...
)::NamedTuple
    effective_config = _effective_transport_config(config, kwargs)

    eta = shear_viscosity(
        quark_params,
        thermo_params,
        effective_config;
        tau=tau,
        degeneracy=degeneracy,
        provider=provider,
    )

    sigma = electric_conductivity(
        quark_params,
        thermo_params,
        effective_config;
        tau=tau,
        charges=charges,
        degeneracy=degeneracy,
        provider=provider,
    )

    zeta = bulk_coeffs === nothing ? NaN : bulk_viscosity_isentropic(
        quark_params,
        thermo_params,
        effective_config;
        tau=tau,
        bulk_coeffs_isentropic=bulk_coeffs,
        provider=provider,
    )

    return (eta=eta, zeta=zeta, sigma=sigma)
end

"""Convenience overload: pass `config` as a positional argument, with optional explicit overrides."""
function transport_coefficients(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    config::TransportIntegrationConfig;
    tau::NamedTuple,
    bulk_coeffs::Union{Nothing,NamedTuple}=nothing,
    charges::NamedTuple=default_charges(),
    degeneracy::Float64=degeneracy_default(),
    provider=DEFAULT_TRANSPORT_PROVIDER,
    kwargs...
)::NamedTuple
    return transport_coefficients(
        quark_params,
        thermo_params;
        tau=tau,
        bulk_coeffs=bulk_coeffs,
        config=config,
        charges=charges,
        degeneracy=degeneracy,
        provider=provider,
        kwargs...
    )
end

end # module
