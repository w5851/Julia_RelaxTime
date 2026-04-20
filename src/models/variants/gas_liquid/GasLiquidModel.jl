"""GasLiquidModel

Gas-Liquid (RMF/Walecka) 最小可运行模型壳。

说明：
- 当前版本先打通 Models 统一接口与 workflow 组装，数值细节在后续 PR 增强。
"""

using StaticArrays

const _GASLIQUID_EQUATION_SET_PATH = normpath(joinpath(@__DIR__, "core", "EquationSet.jl"))
if !isdefined(@__MODULE__, :GasLiquidEquationSet)
    include(_GASLIQUID_EQUATION_SET_PATH)
end

const _GASLIQUID_THERMO_PATH = normpath(joinpath(@__DIR__, "core", "Thermodynamics.jl"))
if !isdefined(@__MODULE__, :GasLiquidThermodynamics)
    include(_GASLIQUID_THERMO_PATH)
end

using .GasLiquidEquationSet: GasLiquidCoreParams, GasLiquidState, solve_equilibrium, mu_baryon
using .GasLiquidThermodynamics: pressure_density_entropy_energy, omega_components as gasliquid_omega_components

export GasLiquidModel

struct GasLiquidModel <: AbstractQCDModel
    params::GasLiquidCoreParams
end

GasLiquidModel(; kwargs...) = GasLiquidModel(GasLiquidCoreParams(; kwargs...))

@inline gap_state_dim(::GasLiquidModel) = 4

function solve_gap(model::GasLiquidModel, T_fm, mu_vec; kwargs...)
    st = solve_equilibrium(T_fm, mu_vec; params=model.params, kwargs...)
    # 统一映射到 MeanFieldState(5) 语义：
    # phi[1]=sigma, phi[2]=delta, phi[3]=0; Phi/PhiBar=1（非 Polyakov 模型）。
    phi = SVector{3, Float64}(st.sigma, st.delta, 0.0)
    return MeanFieldState(phi; Phi=1.0, PhiBar=1.0)
end

@inline function _to_gasliquid_state(x_state, mu_vec)
    st = x_state isa MeanFieldState ? x_state : MeanFieldState(x_state)
    muB = mu_baryon(mu_vec)
    Tprom = promote_type(typeof(st.phi[1]), typeof(st.phi[2]), typeof(muB))
    return GasLiquidState(Tprom(st.phi[1]), Tprom(st.phi[2]), Tprom(muB), Tprom(muB))
end

function calculate_mass_vec(model::GasLiquidModel, φ; kwargs...)
    gl = GasLiquidState(float(φ[1]), float(φ[2]), 0.0, 0.0)
    return GasLiquidThermodynamics.effective_masses(gl, model.params)
end

@inline function calculate_chiral(model::GasLiquidModel, φ; kwargs...)
    _ = kwargs
    sigma = float(φ[1])
    delta = float(φ[2])
    p = model.params
    x = p.g_sigma * sigma
    us = (1 / 3) * p.b * p.m_nucleon_inv_fm * x^3 + (1 / 4) * p.c * x^4
    return us + 0.5 * p.m_delta_inv_fm^2 * delta^2
end

@inline function polyakov_potential(::GasLiquidModel, Φ, Φbar, T; kwargs...)
    _ = (Φ, Φbar, T, kwargs)
    return 0.0
end

@inline function vacuum_contribution(model::GasLiquidModel, masses; kwargs...)
    mu_vec = get(kwargs, :mu_vec, 0.0)
    p_num = Int(get(kwargs, :p_num, 96))
    gl = _state_from_masses_and_mu(model, masses, mu_vec)
    return gasliquid_omega_components(gl, 0.0, model.params; p_num=p_num).vac
end

@inline function _state_from_masses_and_mu(model::GasLiquidModel, masses, mu_vec)
    p = model.params
    mp = float(masses[1])
    mn = float(masses[2])
    muB = mu_baryon(mu_vec)

    sigma = (p.m_nucleon_inv_fm - 0.5 * (mp + mn)) / (p.g_sigma + eps(Float64))
    delta = if abs(p.g_delta) <= eps(Float64)
        0.0
    else
        (mn - mp) / (2 * p.g_delta)
    end

    Tprom = promote_type(typeof(sigma), typeof(delta), typeof(muB))
    return GasLiquidState(Tprom(sigma), Tprom(delta), Tprom(muB), Tprom(muB))
end

@inline function _chi_from_state(model::GasLiquidModel, st::GasLiquidState)
    p = model.params
    x = p.g_sigma * st.sigma
    us = (1 / 3) * p.b * p.m_nucleon_inv_fm * x^3 + (1 / 4) * p.c * x^4
    return us + 0.5 * p.m_delta_inv_fm^2 * st.delta^2
end

@inline function thermal_contribution(model::GasLiquidModel, masses, Φ, Φbar, mu_vec, T; kwargs...)
    _ = (Φ, Φbar)
    p_num = Int(get(kwargs, :p_num, 96))
    gl = _state_from_masses_and_mu(model, masses, mu_vec)
    return gasliquid_omega_components(gl, T, model.params; p_num=p_num).therm
end

function number_densities(model::GasLiquidModel, x_state, T, mu_vec; kwargs...)
    p_num = Int(get(kwargs, :p_num, 96))
    gl = _to_gasliquid_state(x_state, mu_vec)
    thermo = pressure_density_entropy_energy(gl, T, model.params; p_num=p_num)
    # 映射为 quark-like 3 分量契约，保持与 TransportWorkflow 兼容。
    Tm = typeof(thermo.rho)
    q = SVector{3, Tm}(thermo.rho / 3, thermo.rho / 3, thermo.rho / 3)
    aq = SVector{3, Tm}(zero(Tm), zero(Tm), zero(Tm))
    return (quark=q, antiquark=aq)
end
