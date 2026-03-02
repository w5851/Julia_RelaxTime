"""
    Core

PNJL 模型核心计算模块，包含积分和热力学量计算。

## 子模块
- `Integrals`: 积分计算（真空项、热项）
- `Thermodynamics`: 热力学量计算（Ω, P, s, ε, ρ, M）
"""
module Core

using StaticArrays

const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

const _THERMO_FACADE_PATH = normpath(joinpath(@__DIR__, "ThermoFacade.jl"))
const ThermoFacade = IncludeOnce.include_once!(Main, :ThermoFacade, _THERMO_FACADE_PATH)

const _PNJL_CORE_PATH = normpath(joinpath(@__DIR__, "..", "PNJLCore.jl"))
if !isdefined(@__MODULE__, :PNJLCore)
    include(_PNJL_CORE_PATH)
end

const _PNJL_PARAMS_REF = Ref{Any}(nothing)

@inline function _pnjl_params()
    p = _PNJL_PARAMS_REF[]
    if p === nothing
        p = PNJLCore.pnjl_params()
        _PNJL_PARAMS_REF[] = p
    end
    return p
end

@inline calculate_mass_vec(φ) = PNJLCore.calculate_mass_vec(_pnjl_params(), φ)
@inline calculate_mass_vec(x_state::AbstractVector) =
    PNJLCore.calculate_mass_vec(_pnjl_params(), SVector{3}(x_state[1], x_state[2], x_state[3]))

@inline calculate_chiral(φ) = PNJLCore.chiral_potential(_pnjl_params(), φ)
@inline calculate_U(T_fm::Real, Φ::Real, Φbar::Real) = PNJLCore.polyakov_potential(_pnjl_params(), Φ, Φbar, T_fm)

@inline function _safe_log(x)
    min_x = one(x) * 1e-16
    x <= zero(x) && return log(min_x)
    return x < min_x ? log(min_x) : log(x)
end

@inline function calculate_U_derivative_T(T_fm::Real, Φ::Real, Φbar::Real)
    p = _pnjl_params()
    T_ratio = p.T0_inv_fm / T_fm
    Ta = p.a0 + p.a1 * T_ratio + p.a2 * T_ratio^2
    Tb = p.b3 * T_ratio^3
    value = 1 - 6 * Φbar * Φ + 4 * (Φbar^3 + Φ^3) - 3 * (Φbar * Φ)^2
    log_value = _safe_log(value)

    U_over_T4 = -0.5 * Ta * Φbar * Φ + Tb * log_value
    dTa_dT = -p.a1 * p.T0_inv_fm / T_fm^2 - 2 * p.a2 * p.T0_inv_fm^2 / T_fm^3
    dTb_dT = -3 * p.b3 * p.T0_inv_fm^3 / T_fm^4

    return 4 * T_fm^3 * U_over_T4 + T_fm^4 * (dTa_dT * (-0.5 * Φbar * Φ) + dTb_dT * log_value)
end

@inline function calculate_pressure(x_state, mu_vec, T_fm, thermal_nodes, xi)
    return ThermoFacade.calculate_pressure_backend(
        x_state,
        mu_vec,
        T_fm;
        model_kind=:PNJL,
        thermal_nodes=thermal_nodes,
        xi=xi,
    )
end

@inline function calculate_omega(x_state, mu_vec, T_fm, thermal_nodes, xi)
    return ThermoFacade.calculate_omega_backend(
        x_state,
        mu_vec,
        T_fm;
        model_kind=:PNJL,
        thermal_nodes=thermal_nodes,
        xi=xi,
    )
end

@inline function calculate_rho(x_state, mu_vec, T_fm, thermal_nodes, xi)
    return ThermoFacade.calculate_rho_backend(
        x_state,
        mu_vec,
        T_fm;
        model_kind=:PNJL,
        thermal_nodes=thermal_nodes,
        xi=xi,
    )
end

@inline function calculate_thermo(x_state, mu_vec, T_fm, thermal_nodes, xi)
    return ThermoFacade.calculate_thermo_backend(
        x_state,
        mu_vec,
        T_fm;
        model_kind=:PNJL,
        thermal_nodes=thermal_nodes,
        xi=xi,
    )
end

@inline function calculate_number_densities(x_state, mu_vec, T_fm, thermal_nodes, xi)
    return ThermoFacade.calculate_number_densities_backend(
        x_state,
        mu_vec,
        T_fm;
        model_kind=:PNJL,
        thermal_nodes=thermal_nodes,
        xi=xi,
    )
end

const Thermodynamics = @__MODULE__
const ρ0 = ThermoFacade.rho0()

include("Integrals.jl")
include("MagneticIntegrals.jl")
include("MagneticThermodynamics.jl")

using .Integrals
using .MagneticIntegrals
using .MagneticThermodynamics

# 重新导出常用函数
export cached_nodes, vacuum_integral, calculate_energy_sum, calculate_log_sum
export DEFAULT_THETA_COUNT, DEFAULT_MOMENTUM_COUNT
export safe_log, calculate_energy_isotropic, calculate_energy_anisotropic

export calculate_mass_vec, calculate_chiral, calculate_U
export calculate_pressure, calculate_omega
export calculate_rho, calculate_thermo, calculate_number_densities
export ρ0

export QUARK_CHARGE_ABS
export alpha_n, energy_landau, smooth_cutoff, resolve_nmax_from_cutoff
export omega0_flavor_landau, omegat_flavor_landau, density_flavor_landau
export MagneticIMCParams, default_imc_params, coupling_GB
export MagneticConfig, default_magnetic_config
export calculate_magnetic_omega_components, calculate_magnetic_omega
export calculate_magnetic_pressure, calculate_magnetic_rho
export calculate_magnetic_number_densities
export magnetic_nmax_convergence_report

end # module Core

