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

using .GasLiquidEquationSet: GasLiquidCoreParams, GasLiquidState, solve_equilibrium
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
    muB = mu_vec isa Real ? mu_vec : (mu_vec[1] + mu_vec[2] + mu_vec[3]) / 3
    Tprom = promote_type(typeof(st.phi[1]), typeof(st.phi[2]), typeof(muB))
    return GasLiquidState(Tprom(st.phi[1]), Tprom(st.phi[2]), Tprom(muB), Tprom(muB))
end

function calculate_mass_vec(model::GasLiquidModel, φ; kwargs...)
    gl = GasLiquidState(float(φ[1]), float(φ[2]), 0.0, 0.0)
    return GasLiquidThermodynamics.effective_masses(gl, model.params)
end

@inline function calculate_chiral(::GasLiquidModel, φ; kwargs...)
    return 0.05 * (float(φ[1])^2 + float(φ[2])^2)
end

@inline function polyakov_potential(::GasLiquidModel, Φ, Φbar, T; kwargs...)
    _ = (Φ, Φbar, T, kwargs)
    return 0.0
end

@inline function vacuum_contribution(model::GasLiquidModel, masses; kwargs...)
    _ = (kwargs,)
    return -0.03 * sum(float.(masses))
end

@inline function thermal_contribution(model::GasLiquidModel, masses, Φ, Φbar, mu_vec, T; kwargs...)
    _ = (model, masses, Φ, Φbar, kwargs)
    muB = mu_vec isa Real ? float(mu_vec) : float((mu_vec[1] + mu_vec[2] + mu_vec[3]) / 3)
    return -0.02 * float(T) * (1 + abs(muB))
end

function number_densities(model::GasLiquidModel, x_state, T, mu_vec; kwargs...)
    _ = (kwargs,)
    gl = _to_gasliquid_state(x_state, mu_vec)
    thermo = pressure_density_entropy_energy(gl, T, model.params)
    # 映射为 quark-like 3 分量契约，保持与 TransportWorkflow 兼容。
    Tm = typeof(thermo.rho)
    q = SVector{3, Tm}(thermo.rho, thermo.rho, thermo.rho)
    aq = SVector{3, Tm}(zero(Tm), zero(Tm), zero(Tm))
    return (quark=q, antiquark=aq)
end

function omega_components(model::GasLiquidModel, x_state, T, mu_vec; kwargs...)
    _ = (kwargs,)
    gl = _to_gasliquid_state(x_state, mu_vec)
    return gasliquid_omega_components(gl, T, model.params)
end
