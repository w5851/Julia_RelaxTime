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

# Dependencies loaded by RelaxTime.jl entry point

using Base.MathConstants: π

using Main.Constants_PNJL: N_color, α
using ..GaussLegendre: gauleg, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS, DEFAULT_COSΘ_NODES, DEFAULT_COSΘ_WEIGHTS
using ..PhaseSpaceSampling: p_nodes_weights, cos_nodes_weights, integrate_p, integrate_p_cos
using Main.PNJLQuarkDistributions: quark_distribution, antiquark_distribution
using Main.PNJLQuarkDistributions_Aniso: quark_distribution_aniso, antiquark_distribution_aniso

using Main.ParameterTypes: QuarkParams, ThermoParams
using Main.ParameterAdapters: normalize_quark_input, normalize_thermo_input
using Main.ValidationUtils: require_finite, require_positive_finite,
    require_nonnegative_finite, require_positive_integer, validate_grid_weight_pair

export shear_viscosity, bulk_viscosity, bulk_viscosity_isentropic, electric_conductivity, transport_coefficients
export diffusion_coefficient, diffusion_matrix
export kappa_BB, kappa_BQ, kappa_BS, kappa_QQ, kappa_QS, kappa_SS, lambda_from_kappa_BB
export conserved_charge_value, conserved_charge_densities, enthalpy_density, rho_mass_from_densities
export lorenz_number, lorentz_legacy, viscous_conductive_coupling_ratio, prandtl_number, bulk_to_shear_viscosity_ratio

export TransportIntegrationConfig, QuarkParams, ThermoParams, TransportPhysicsConfig, TransportRequest
export default_transport_provider

const TWO_PI = 2.0 * π
const ENERGY_FLOOR = sqrt(eps(Float64))

const _SPECIES_ALL = (:u, :d, :s, :ubar, :dbar, :sbar)
const _FLAVORS = (:u, :d, :s)
const _CONSERVED_CHARGES = (:B, :Q, :S)

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
    require_positive_integer("TransportIntegrationConfig.p_nodes", p_nodes)
    require_positive_integer("TransportIntegrationConfig.cos_nodes", cos_nodes)
    require_positive_finite("TransportIntegrationConfig.p_max", p_max)
    validate_grid_weight_pair("TransportIntegrationConfig", "p_grid", p_grid, "p_w", p_w)
    validate_grid_weight_pair("TransportIntegrationConfig", "cos_grid", cos_grid, "cos_w", cos_w)
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

@inline _qp_view(q::QuarkParams) = normalize_quark_input(q)
@inline _tp_view(t::ThermoParams) = normalize_thermo_input(t)

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
    densities::Union{Nothing,NamedTuple}=nothing,
    pressure::Union{Nothing,Real}=nothing,
    energy::Union{Nothing,Real}=nothing,
    entropy::Union{Nothing,Real}=nothing,
    c_p::Union{Nothing,Real}=nothing,
    rho_mass::Union{Nothing,Real}=nothing,
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
        densities=densities,
        pressure=pressure,
        energy=energy,
        entropy=entropy,
        c_p=c_p,
        rho_mass=rho_mass,
        charges=charges,
        degeneracy=degeneracy,
        config=integration,
        provider=provider,
        kwargs...
    )
end

@inline function fermi_factor(f::Float64)
    return f * (1.0 - f)
end

@inline function conserved_charge_value(species::Symbol, charge::Symbol)::Float64
    if charge === :B
        return _is_quark_species(species) ? (1.0 / 3.0) : (-1.0 / 3.0)
    elseif charge === :Q
        if species === :u
            return 2.0 / 3.0
        elseif species === :d || species === :s
            return -1.0 / 3.0
        elseif species === :ubar
            return -2.0 / 3.0
        elseif species === :dbar || species === :sbar
            return 1.0 / 3.0
        end
    elseif charge === :S
        if species === :s
            return -1.0
        elseif species === :sbar
            return 1.0
        elseif species === :u || species === :d || species === :ubar || species === :dbar
            return 0.0
        end
    end
    error("Unknown conserved charge or species: charge=$charge species=$species")
end

function conserved_charge_densities(densities::NamedTuple)::NamedTuple
    for sp in _SPECIES_ALL
        hasproperty(densities, sp) || error("densities is missing :$sp")
        val = getproperty(densities, sp)
        isfinite(val) || error("densities.:$sp must be finite")
    end

    Δu = Float64(densities.u - densities.ubar)
    Δd = Float64(densities.d - densities.dbar)
    Δs = Float64(densities.s - densities.sbar)

    return (
        B=(Δu + Δd + Δs) / 3.0,
        Q=(2.0 * Δu - Δd - Δs) / 3.0,
        S=-Δs,
    )
end

@inline function enthalpy_density(pressure::Real, energy::Real)::Float64
    p = Float64(pressure)
    e = Float64(energy)
    isfinite(p) || error("pressure must be finite")
    isfinite(e) || error("energy must be finite")
    return p + e
end

function rho_mass_from_densities(masses::NamedTuple, densities::NamedTuple)::Float64
    for flavor in _FLAVORS
        hasproperty(masses, flavor) || error("masses is missing :$flavor")
        mass = Float64(getproperty(masses, flavor))
        isfinite(mass) || error("masses.:$flavor must be finite")
        mass >= 0.0 || error("masses.:$flavor must be >= 0")
    end
    for sp in _SPECIES_ALL
        hasproperty(densities, sp) || error("densities is missing :$sp")
        val = Float64(getproperty(densities, sp))
        isfinite(val) || error("densities.:$sp must be finite")
    end

    return (Float64(densities.u) + Float64(densities.ubar)) * Float64(masses.u) +
           (Float64(densities.d) + Float64(densities.dbar)) * Float64(masses.d) +
           (Float64(densities.s) + Float64(densities.sbar)) * Float64(masses.s)
end

@inline function _safe_ratio(numerator::Real, denominator::Real)::Float64
    num = Float64(numerator)
    den = Float64(denominator)
    isfinite(num) || return NaN
    isfinite(den) || return NaN
    abs(den) <= sqrt(eps(Float64)) && return NaN
    return num / den
end

@inline function _sigma_over_temperature(sigma::Real, T::Real)::Float64
    temp = Float64(T)
    isfinite(temp) && temp > 0.0 || error("T must be finite and > 0")
    return _safe_ratio(sigma, temp)
end

@inline function lorenz_number(lambda_value::Real, sigma::Real, T::Real)::Float64
    temp = Float64(T)
    isfinite(temp) && temp > 0.0 || error("T must be finite and > 0")
    return _safe_ratio(lambda_value, Float64(sigma) * temp)
end

@inline function lorentz_legacy(lambda_value::Real, sigma::Real, T::Real)::Float64
    return _safe_ratio(lambda_value, _sigma_over_temperature(sigma, T))
end

@inline function viscous_conductive_coupling_ratio(eta::Real, entropy_density::Real, sigma::Real, T::Real)::Float64
    eta_over_s = _safe_ratio(eta, entropy_density)
    return _safe_ratio(eta_over_s, _sigma_over_temperature(sigma, T))
end

@inline function prandtl_number(eta::Real, c_p::Real, lambda_value::Real, rho_mass::Real)::Float64
    return _safe_ratio(Float64(eta) * Float64(c_p), Float64(lambda_value) * Float64(rho_mass))
end

@inline function bulk_to_shear_viscosity_ratio(zeta::Real, eta::Real)::Float64
    return _safe_ratio(zeta, eta)
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
)::NTuple{4,Float64}
    T = thermo_params.T
    Φ = thermo_params.Φ
    Φbar = thermo_params.Φbar

    m = _mass_for_species(provider, caps, sp, quark_params, thermo_params)
    μ = _mu_for_species(provider, caps, sp, quark_params, thermo_params)
    E = _energy_for_kernel(provider, caps, p, m, ξ, c)
    f = distribution_for_species_from_E(provider, caps, sp, E, p, m, μ, T, Φ, Φbar, ξ, c)
    ff = fermi_factor(f)
    τ = tau_for_species(sp, tau)
    return (E, f, ff, τ)
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
    require_positive_finite("thermo_params.T", T)
    require_finite("thermo_params.Φ", thermo_params.Φ)
    require_finite("thermo_params.Φbar", thermo_params.Φbar)

    for sp in _SPECIES_ALL
        τ = tau_for_species(sp, tau)
        require_nonnegative_finite("tau.:$sp", τ)
    end

    check_mass = provider === nothing || !hasproperty(provider, :mass_for_species)
    check_mu = provider === nothing || !hasproperty(provider, :mu_for_species)

    for fl in _FLAVORS
        m = getproperty(quark_params.m, fl)
        μ = getproperty(quark_params.μ, fl)
        if check_mass
            require_nonnegative_finite("quark_params.m.$fl", m)
        end
        if check_mu
            require_finite("quark_params.μ.$fl", μ)
        end
    end

    require_positive_integer("integration p_nodes", config.p_nodes)
    require_positive_integer("integration cos_nodes", config.cos_nodes)
    require_positive_finite("integration p_max", config.p_max)

    if charges !== nothing
        require_finite("charges.u", charges.u)
        require_finite("charges.d", charges.d)
        require_finite("charges.s", charges.s)
    end

    if bulk_coeffs_isentropic !== nothing
        require_finite("bulk_coeffs_isentropic.v_n_sq", bulk_coeffs_isentropic.v_n_sq)
        require_finite("bulk_coeffs_isentropic.dμB_dT_sigma", bulk_coeffs_isentropic.dμB_dT_sigma)

        length(bulk_coeffs_isentropic.masses) == 3 || throw(ArgumentError("bulk_coeffs_isentropic.masses must have length 3"))
        length(bulk_coeffs_isentropic.dM_dT) == 3 || throw(ArgumentError("bulk_coeffs_isentropic.dM_dT must have length 3"))
        length(bulk_coeffs_isentropic.dM_dμB) == 3 || throw(ArgumentError("bulk_coeffs_isentropic.dM_dμB must have length 3"))

        all(isfinite, bulk_coeffs_isentropic.masses) || throw(ArgumentError("bulk_coeffs_isentropic.masses must be finite"))
        if !all(m -> m >= 0.0, bulk_coeffs_isentropic.masses)
            @warn "bulk_coeffs_isentropic.masses contains negative value(s) — may indicate unphysical gap solution; proceeding with abs(m)" masses=bulk_coeffs_isentropic.masses
        end
        all(isfinite, bulk_coeffs_isentropic.dM_dT) || throw(ArgumentError("bulk_coeffs_isentropic.dM_dT must be finite"))
        all(isfinite, bulk_coeffs_isentropic.dM_dμB) || throw(ArgumentError("bulk_coeffs_isentropic.dM_dμB must be finite"))
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
                E, f, ff, τ = _species_transport_state(provider, caps, sp, p, 0.0, 0.0, quark_params, thermo_params, tau)
                inner += term_for_species(p, sp, E, f, ff, τ)
            end
            return inner
        end
        return pref_measure_iso * acc
    end

    nodes_cos, weights_cos = _cos_quadrature(config)
    acc = integrate_p_cos(nodes_p, weights_p, nodes_cos, weights_cos) do p, c
        inner = 0.0
        for sp in _SPECIES_ALL
            E, f, ff, τ = _species_transport_state(provider, caps, sp, p, c, ξ, quark_params, thermo_params, tau)
            inner += term_for_species(p, sp, E, f, ff, τ)
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

    integral = _integrate_species_sum(quark_params, thermo_params, tau, provider, effective_config) do p, sp, E, f, ff, τ
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

    integral = _integrate_species_sum(quark_params, thermo_params, tau, provider, effective_config) do p, sp, E, f, ff, τ
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
    masses = bulk_coeffs_isentropic.masses  # 使用模型给出的质量（可能为负，物理上保留）
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

@inline function _validate_conserved_charge_symbol(charge::Symbol)::Symbol
    charge in _CONSERVED_CHARGES || error("Unknown conserved charge: $charge. Allowed charges: $(_CONSERVED_CHARGES)")
    return charge
end

function _validate_diffusion_background(densities::NamedTuple, pressure::Real, energy::Real)
    charge_densities = conserved_charge_densities(densities)
    h = enthalpy_density(pressure, energy)
    h > 0.0 || error("enthalpy density must be > 0")
    return charge_densities, h
end

@inline function _kappa_projection(species::Symbol, charge::Symbol, charge_densities::NamedTuple, enthalpy::Float64, E::Float64)::Float64
    q = conserved_charge_value(species, charge)
    n = getproperty(charge_densities, charge)
    return q - (n / enthalpy) * E
end

@inline function _kappa_integrand_weight(
    species::Symbol,
    charge_left::Symbol,
    charge_right::Symbol,
    charge_densities::NamedTuple,
    enthalpy::Float64,
    E::Float64,
)::Float64
    left = _kappa_projection(species, charge_left, charge_densities, enthalpy, E)
    right = _kappa_projection(species, charge_right, charge_densities, enthalpy, E)
    return left * right
end

function diffusion_coefficient(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    tau::NamedTuple,
    charge_left::Symbol,
    charge_right::Symbol,
    densities::NamedTuple,
    pressure::Real,
    energy::Real,
    degeneracy::Float64=degeneracy_default(),
    provider=DEFAULT_TRANSPORT_PROVIDER,
    config::Union{Nothing,TransportIntegrationConfig}=nothing,
    kwargs...
)::Float64
    T = thermo_params.T
    charge_left = _validate_conserved_charge_symbol(charge_left)
    charge_right = _validate_conserved_charge_symbol(charge_right)

    effective_config = _effective_transport_config(config, kwargs)
    _validate_transport_inputs(quark_params, thermo_params, tau, effective_config; provider=provider)
    charge_densities, h = _validate_diffusion_background(densities, pressure, energy)

    integral = _integrate_species_sum(quark_params, thermo_params, tau, provider, effective_config) do p, sp, E, f, ff, τ
        p2 = p * p
        p4 = p2 * p2
        weight = _kappa_integrand_weight(sp, charge_left, charge_right, charge_densities, h, E)
        return p4 / (E * E) * (degeneracy * τ * f * weight)
    end

    return integral / (3.0 * T * T)
end

function diffusion_coefficient(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    config::TransportIntegrationConfig;
    tau::NamedTuple,
    charge_left::Symbol,
    charge_right::Symbol,
    densities::NamedTuple,
    pressure::Real,
    energy::Real,
    degeneracy::Float64=degeneracy_default(),
    provider=DEFAULT_TRANSPORT_PROVIDER,
    kwargs...
)::Float64
    return diffusion_coefficient(
        quark_params,
        thermo_params;
        tau=tau,
        charge_left=charge_left,
        charge_right=charge_right,
        densities=densities,
        pressure=pressure,
        energy=energy,
        degeneracy=degeneracy,
        provider=provider,
        config=config,
        kwargs...
    )
end

function diffusion_coefficient(
    quark::QuarkParams,
    thermo::ThermoParams;
    kwargs...
)::Float64
    return diffusion_coefficient(_qp_view(quark), _tp_view(thermo); kwargs...)
end

function diffusion_coefficient(
    quark::QuarkParams,
    thermo::ThermoParams,
    config::TransportIntegrationConfig;
    kwargs...
)::Float64
    return diffusion_coefficient(_qp_view(quark), _tp_view(thermo), config; kwargs...)
end

function diffusion_coefficient(
    req::TransportRequest;
    tau::NamedTuple=req.tau,
    charge_left::Symbol,
    charge_right::Symbol,
    densities::NamedTuple,
    pressure::Real,
    energy::Real,
    degeneracy::Float64=req.physics.degeneracy,
    integration::TransportIntegrationConfig=req.integration,
    provider=DEFAULT_TRANSPORT_PROVIDER,
    kwargs...
)::Float64
    quark_params = (m=req.quark.m, μ=req.quark.μ)
    thermo_params = (T=req.thermo.T, Φ=req.thermo.Φ, Φbar=req.thermo.Φbar, ξ=req.thermo.ξ)
    return diffusion_coefficient(
        quark_params,
        thermo_params;
        tau=tau,
        charge_left=charge_left,
        charge_right=charge_right,
        densities=densities,
        pressure=pressure,
        energy=energy,
        degeneracy=degeneracy,
        config=integration,
        provider=provider,
        kwargs...
    )
end

@inline function kappa_BB(quark_params::NamedTuple, thermo_params::NamedTuple; kwargs...)::Float64
    return diffusion_coefficient(quark_params, thermo_params; charge_left=:B, charge_right=:B, kwargs...)
end

@inline function kappa_BB(quark::QuarkParams, thermo::ThermoParams; kwargs...)::Float64
    return kappa_BB(_qp_view(quark), _tp_view(thermo); kwargs...)
end

@inline function kappa_BQ(quark_params::NamedTuple, thermo_params::NamedTuple; kwargs...)::Float64
    return diffusion_coefficient(quark_params, thermo_params; charge_left=:B, charge_right=:Q, kwargs...)
end

@inline function kappa_BQ(quark::QuarkParams, thermo::ThermoParams; kwargs...)::Float64
    return kappa_BQ(_qp_view(quark), _tp_view(thermo); kwargs...)
end

@inline function kappa_BS(quark_params::NamedTuple, thermo_params::NamedTuple; kwargs...)::Float64
    return diffusion_coefficient(quark_params, thermo_params; charge_left=:B, charge_right=:S, kwargs...)
end

@inline function kappa_BS(quark::QuarkParams, thermo::ThermoParams; kwargs...)::Float64
    return kappa_BS(_qp_view(quark), _tp_view(thermo); kwargs...)
end

@inline function kappa_QQ(quark_params::NamedTuple, thermo_params::NamedTuple; kwargs...)::Float64
    return diffusion_coefficient(quark_params, thermo_params; charge_left=:Q, charge_right=:Q, kwargs...)
end

@inline function kappa_QQ(quark::QuarkParams, thermo::ThermoParams; kwargs...)::Float64
    return kappa_QQ(_qp_view(quark), _tp_view(thermo); kwargs...)
end

@inline function kappa_QS(quark_params::NamedTuple, thermo_params::NamedTuple; kwargs...)::Float64
    return diffusion_coefficient(quark_params, thermo_params; charge_left=:Q, charge_right=:S, kwargs...)
end

@inline function kappa_QS(quark::QuarkParams, thermo::ThermoParams; kwargs...)::Float64
    return kappa_QS(_qp_view(quark), _tp_view(thermo); kwargs...)
end

@inline function kappa_SS(quark_params::NamedTuple, thermo_params::NamedTuple; kwargs...)::Float64
    return diffusion_coefficient(quark_params, thermo_params; charge_left=:S, charge_right=:S, kwargs...)
end

@inline function kappa_SS(quark::QuarkParams, thermo::ThermoParams; kwargs...)::Float64
    return kappa_SS(_qp_view(quark), _tp_view(thermo); kwargs...)
end

function diffusion_matrix(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    tau::NamedTuple,
    densities::NamedTuple,
    pressure::Real,
    energy::Real,
    degeneracy::Float64=degeneracy_default(),
    provider=DEFAULT_TRANSPORT_PROVIDER,
    config::Union{Nothing,TransportIntegrationConfig}=nothing,
    kwargs...
)::NamedTuple
    κBB = kappa_BB(
        quark_params,
        thermo_params;
        tau=tau,
        densities=densities,
        pressure=pressure,
        energy=energy,
        degeneracy=degeneracy,
        provider=provider,
        config=config,
        kwargs...
    )
    κBQ = kappa_BQ(
        quark_params,
        thermo_params;
        tau=tau,
        densities=densities,
        pressure=pressure,
        energy=energy,
        degeneracy=degeneracy,
        provider=provider,
        config=config,
        kwargs...
    )
    κBS = kappa_BS(
        quark_params,
        thermo_params;
        tau=tau,
        densities=densities,
        pressure=pressure,
        energy=energy,
        degeneracy=degeneracy,
        provider=provider,
        config=config,
        kwargs...
    )
    κQQ = kappa_QQ(
        quark_params,
        thermo_params;
        tau=tau,
        densities=densities,
        pressure=pressure,
        energy=energy,
        degeneracy=degeneracy,
        provider=provider,
        config=config,
        kwargs...
    )
    κQS = kappa_QS(
        quark_params,
        thermo_params;
        tau=tau,
        densities=densities,
        pressure=pressure,
        energy=energy,
        degeneracy=degeneracy,
        provider=provider,
        config=config,
        kwargs...
    )
    κSS = kappa_SS(
        quark_params,
        thermo_params;
        tau=tau,
        densities=densities,
        pressure=pressure,
        energy=energy,
        degeneracy=degeneracy,
        provider=provider,
        config=config,
        kwargs...
    )

    return (
        charges=_CONSERVED_CHARGES,
        matrix=[κBB κBQ κBS; κBQ κQQ κQS; κBS κQS κSS],
        BB=κBB,
        BQ=κBQ,
        BS=κBS,
        QB=κBQ,
        QQ=κQQ,
        QS=κQS,
        SB=κBS,
        SQ=κQS,
        SS=κSS,
    )
end

function diffusion_matrix(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    config::TransportIntegrationConfig;
    tau::NamedTuple,
    densities::NamedTuple,
    pressure::Real,
    energy::Real,
    degeneracy::Float64=degeneracy_default(),
    provider=DEFAULT_TRANSPORT_PROVIDER,
    kwargs...
)::NamedTuple
    return diffusion_matrix(
        quark_params,
        thermo_params;
        tau=tau,
        densities=densities,
        pressure=pressure,
        energy=energy,
        degeneracy=degeneracy,
        provider=provider,
        config=config,
        kwargs...
    )
end

function diffusion_matrix(
    quark::QuarkParams,
    thermo::ThermoParams;
    kwargs...
)::NamedTuple
    return diffusion_matrix(_qp_view(quark), _tp_view(thermo); kwargs...)
end

function diffusion_matrix(
    quark::QuarkParams,
    thermo::ThermoParams,
    config::TransportIntegrationConfig;
    kwargs...
)::NamedTuple
    return diffusion_matrix(_qp_view(quark), _tp_view(thermo), config; kwargs...)
end

function diffusion_matrix(
    req::TransportRequest;
    tau::NamedTuple=req.tau,
    densities::NamedTuple,
    pressure::Real,
    energy::Real,
    degeneracy::Float64=req.physics.degeneracy,
    integration::TransportIntegrationConfig=req.integration,
    provider=DEFAULT_TRANSPORT_PROVIDER,
    kwargs...
)::NamedTuple
    quark_params = (m=req.quark.m, μ=req.quark.μ)
    thermo_params = (T=req.thermo.T, Φ=req.thermo.Φ, Φbar=req.thermo.Φbar, ξ=req.thermo.ξ)
    return diffusion_matrix(
        quark_params,
        thermo_params;
        tau=tau,
        densities=densities,
        pressure=pressure,
        energy=energy,
        degeneracy=degeneracy,
        config=integration,
        provider=provider,
        kwargs...
    )
end

@inline function kappa_BB(req::TransportRequest; kwargs...)::Float64
    return kappa_BB((m=req.quark.m, μ=req.quark.μ), (T=req.thermo.T, Φ=req.thermo.Φ, Φbar=req.thermo.Φbar, ξ=req.thermo.ξ); tau=req.tau, kwargs...)
end

@inline function kappa_BQ(req::TransportRequest; kwargs...)::Float64
    return kappa_BQ((m=req.quark.m, μ=req.quark.μ), (T=req.thermo.T, Φ=req.thermo.Φ, Φbar=req.thermo.Φbar, ξ=req.thermo.ξ); tau=req.tau, kwargs...)
end

@inline function kappa_BS(req::TransportRequest; kwargs...)::Float64
    return kappa_BS((m=req.quark.m, μ=req.quark.μ), (T=req.thermo.T, Φ=req.thermo.Φ, Φbar=req.thermo.Φbar, ξ=req.thermo.ξ); tau=req.tau, kwargs...)
end

@inline function kappa_QQ(req::TransportRequest; kwargs...)::Float64
    return kappa_QQ((m=req.quark.m, μ=req.quark.μ), (T=req.thermo.T, Φ=req.thermo.Φ, Φbar=req.thermo.Φbar, ξ=req.thermo.ξ); tau=req.tau, kwargs...)
end

@inline function kappa_QS(req::TransportRequest; kwargs...)::Float64
    return kappa_QS((m=req.quark.m, μ=req.quark.μ), (T=req.thermo.T, Φ=req.thermo.Φ, Φbar=req.thermo.Φbar, ξ=req.thermo.ξ); tau=req.tau, kwargs...)
end

@inline function kappa_SS(req::TransportRequest; kwargs...)::Float64
    return kappa_SS((m=req.quark.m, μ=req.quark.μ), (T=req.thermo.T, Φ=req.thermo.Φ, Φbar=req.thermo.Φbar, ξ=req.thermo.ξ); tau=req.tau, kwargs...)
end

function lambda_from_kappa_BB(kappa_BB::Real, pressure::Real, energy::Real, n_B::Real, T::Real)::Float64
    κ = Float64(kappa_BB)
    h = enthalpy_density(pressure, energy)
    nB = Float64(n_B)
    temp = Float64(T)
    isfinite(κ) || return NaN
    isfinite(nB) || error("n_B must be finite")
    isfinite(temp) && temp > 0.0 || error("T must be finite and > 0")
    if abs(nB) <= sqrt(eps(Float64))
        return NaN
    end
    return κ * (h / (nB * temp))^2
end

function lambda_from_kappa_BB(quark_params::NamedTuple, thermo_params::NamedTuple; kwargs...)::Float64
    haskey(kwargs, :densities) || error("lambda_from_kappa_BB requires densities")
    haskey(kwargs, :pressure) || error("lambda_from_kappa_BB requires pressure")
    haskey(kwargs, :energy) || error("lambda_from_kappa_BB requires energy")
    densities = kwargs[:densities]
    pressure = kwargs[:pressure]
    energy = kwargs[:energy]
    κBB = kappa_BB(quark_params, thermo_params; kwargs...)
    nB = conserved_charge_densities(densities).B
    return lambda_from_kappa_BB(κBB, pressure, energy, nB, thermo_params.T)
end

function lambda_from_kappa_BB(quark::QuarkParams, thermo::ThermoParams; kwargs...)::Float64
    return lambda_from_kappa_BB(_qp_view(quark), _tp_view(thermo); kwargs...)
end

function lambda_from_kappa_BB(req::TransportRequest; kwargs...)::Float64
    return lambda_from_kappa_BB((m=req.quark.m, μ=req.quark.μ), (T=req.thermo.T, Φ=req.thermo.Φ, Φbar=req.thermo.Φbar, ξ=req.thermo.ξ); tau=req.tau, kwargs...)
end

"""
    transport_coefficients(quark_params, thermo_params; tau, bulk_coeffs, charges, ...)

一次性计算 (η, ζ, σ, κ 对角元/非对角元, λ) 及其派生诊断比值。

- `lorenz_number = lambda / (sigma * T)` 是标准 Lorenz 数。
- `lorentz_legacy = lambda / (sigma / T)` 保留 legacy Fortran 诊断口径。
- `viscous_conductive_coupling_ratio` 当前沿用 legacy 归一化 `(eta / s) / (sigma / T)`。
- `prandtl_number` 需要额外提供 `c_p` 与 `rho_mass`。
- `bulk_to_shear_viscosity_ratio = zeta / eta`。
"""
function transport_coefficients(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    tau::NamedTuple,
    bulk_coeffs::Union{Nothing,NamedTuple}=nothing,
    densities::Union{Nothing,NamedTuple}=nothing,
    pressure::Union{Nothing,Real}=nothing,
    energy::Union{Nothing,Real}=nothing,
    entropy::Union{Nothing,Real}=nothing,
    c_p::Union{Nothing,Real}=nothing,
    rho_mass::Union{Nothing,Real}=nothing,
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

    has_diffusion_background = densities !== nothing || pressure !== nothing || energy !== nothing
    if has_diffusion_background
        densities !== nothing || error("transport_coefficients requires densities when pressure/energy is provided")
        pressure !== nothing || error("transport_coefficients requires pressure when densities/energy is provided")
        energy !== nothing || error("transport_coefficients requires energy when densities/pressure is provided")
    end

    κmat = has_diffusion_background ? diffusion_matrix(
        quark_params,
        thermo_params,
        effective_config;
        tau=tau,
        densities=densities,
        pressure=pressure,
        energy=energy,
        degeneracy=degeneracy,
        provider=provider,
    ) : nothing

    κBB = κmat === nothing ? NaN : κmat.BB
    κBQ = κmat === nothing ? NaN : κmat.BQ
    κBS = κmat === nothing ? NaN : κmat.BS
    κQQ = κmat === nothing ? NaN : κmat.QQ
    κQS = κmat === nothing ? NaN : κmat.QS
    κSS = κmat === nothing ? NaN : κmat.SS

    λ = has_diffusion_background ? lambda_from_kappa_BB(κBB, pressure, energy, conserved_charge_densities(densities).B, thermo_params.T) : NaN
    L = lorenz_number(λ, sigma, thermo_params.T)
    L_legacy = lorentz_legacy(λ, sigma, thermo_params.T)
    ratio_eta_sigma = entropy === nothing ? NaN : viscous_conductive_coupling_ratio(eta, entropy, sigma, thermo_params.T)
    Pr = (c_p === nothing || rho_mass === nothing) ? NaN : prandtl_number(eta, c_p, λ, rho_mass)
    zeta_over_eta = bulk_to_shear_viscosity_ratio(zeta, eta)

    return (
        eta=eta,
        zeta=zeta,
        sigma=sigma,
        kappa_BB=κBB,
        kappa_BQ=κBQ,
        kappa_BS=κBS,
        kappa_QQ=κQQ,
        kappa_QS=κQS,
        kappa_SS=κSS,
        diffusion_matrix=κmat,
        lambda=λ,
        lorenz_number=L,
        lorentz_legacy=L_legacy,
        viscous_conductive_coupling_ratio=ratio_eta_sigma,
        prandtl_number=Pr,
        bulk_to_shear_viscosity_ratio=zeta_over_eta,
    )
end

"""Convenience overload: pass `config` as a positional argument, with optional explicit overrides."""
function transport_coefficients(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    config::TransportIntegrationConfig;
    tau::NamedTuple,
    bulk_coeffs::Union{Nothing,NamedTuple}=nothing,
    densities::Union{Nothing,NamedTuple}=nothing,
    pressure::Union{Nothing,Real}=nothing,
    energy::Union{Nothing,Real}=nothing,
    entropy::Union{Nothing,Real}=nothing,
    c_p::Union{Nothing,Real}=nothing,
    rho_mass::Union{Nothing,Real}=nothing,
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
        densities=densities,
        pressure=pressure,
        energy=energy,
        entropy=entropy,
        c_p=c_p,
        rho_mass=rho_mass,
        config=config,
        charges=charges,
        degeneracy=degeneracy,
        provider=provider,
        kwargs...
    )
end

end # module
