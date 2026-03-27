"""PNJLMagneticModel

models 侧外磁场 PNJL 适配器（最小可用版）。
"""

using StaticArrays

export PNJLMagneticModel

const _PNJL_MAG_CORE_PATH = normpath(joinpath(@__DIR__, "core", "MagneticThermodynamics.jl"))
if !isdefined(Main, :MagneticThermodynamics)
    Base.include(Main, _PNJL_MAG_CORE_PATH)
end

@inline function _magnetic_thermodynamics_module()
    if isdefined(Main, :Models) && isdefined(Main.Models, :magnetic_thermodynamics_module)
        return Main.Models.magnetic_thermodynamics_module()
    end
    return Main.MagneticThermodynamics
end

struct PNJLMagneticModel{MT} <: AbstractPNJLModel
    base::PNJLModel
    magnetic::MT
end

function PNJLMagneticModel(; eB_fm2::Real=0.0, profile::String=get(ENV, "PNJL_PARAM_PROFILE", "default"), physics_profile::String=get(ENV, "PHYSICS_PARAM_PROFILE", "default"), kwargs...)
    base = PNJLModel(; profile=profile, physics_profile=physics_profile)
    conf = _magnetic_thermodynamics_module().default_magnetic_config(eB_fm2=float(eB_fm2))
    return PNJLMagneticModel(base, conf)
end

@inline function calculate_mass_vec(model::PNJLMagneticModel, φ; kwargs...)
    return calculate_mass_vec(model.base, φ; kwargs...)
end

@inline function calculate_chiral(model::PNJLMagneticModel, φ; kwargs...)
    return calculate_chiral(model.base, φ; kwargs...)
end

@inline function vacuum_contribution(model::PNJLMagneticModel, masses; kwargs...)
    return vacuum_contribution(model.base, masses; kwargs...)
end

@inline function thermal_contribution(model::PNJLMagneticModel, masses, Φ, Φbar, mu_vec, T; kwargs...)
    return thermal_contribution(model.base, masses, Φ, Φbar, mu_vec, T; kwargs...)
end

@inline function polyakov_potential(model::PNJLMagneticModel, Φ, Φbar, T; kwargs...)
    return polyakov_potential(model.base, Φ, Φbar, T; kwargs...)
end

@inline function number_densities(model::PNJLMagneticModel, x_state, T, mu_vec; kwargs...)
    return _magnetic_thermodynamics_module().calculate_magnetic_number_densities(x_state, normalize_mu_vec(mu_vec), T, model.magnetic)
end

function solve_gap(model::PNJLMagneticModel, T_fm, mu_vec; kwargs...)
    return solve_gap(model.base, T_fm, mu_vec; kwargs...)
end

function omega_components(model::PNJLMagneticModel, x_state, T, mu_vec; kwargs...)
    st = x_state isa MeanFieldState ? state_vector(x_state) : SVector{5, Float64}(Tuple(x_state))
    μ = normalize_mu_vec(mu_vec)
    comp = _magnetic_thermodynamics_module().calculate_magnetic_omega_components(st, μ, T, model.magnetic)
    return (chi=comp.chi, poly=comp.poly, vac=comp.vac, therm=comp.therm, masses=comp.masses, omega=comp.omega)
end

@inline function omega(model::PNJLMagneticModel, x_state, T, mu_vec; kwargs...)
    return omega_components(model, x_state, T, mu_vec; kwargs...).omega
end
