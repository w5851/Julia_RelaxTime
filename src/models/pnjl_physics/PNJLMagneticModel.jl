"""PNJLMagneticModel

models 侧外磁场 PNJL 适配器（最小可用版）。
"""

using StaticArrays
import ..Models

export PNJLMagneticModel
export MagneticGapCandidate, MagneticGapResult
export magnetic_gap_residual, magnetic_gap_residual_autodiff, solve_magnetic_gap

@inline function _magnetic_thermodynamics_module()
    if isdefined(Models, :MagneticThermodynamics)
        return Models.MagneticThermodynamics
    end
    return Main.MagneticThermodynamics
end

struct PNJLMagneticModel{MT} <: AbstractPNJLModel
    base::PNJLModel
    magnetic::MT
end

include(joinpath(@__DIR__, "MagneticGapSolver.jl"))

@inline function model_capabilities(model::PNJLMagneticModel)
    # The generic density contract requires independent quark and antiquark
    # vectors.  The magnetic route exposes only the dedicated net-density API.
    return ModelCapabilities(
        supports_solve_gap=true,
        supports_model_thermo=true,
        supports_number_densities=false,
    )
end

function PNJLMagneticModel(
    ;
    eB_fm2::Real=0.0,
    profile::String=get(ENV, "PNJL_PARAM_PROFILE", "default"),
    physics_profile::String=get(ENV, "PHYSICS_PARAM_PROFILE", "default"),
    magnetic_profile::String=get(ENV, "PNJL_MAGNETIC_PROFILE", "magnetic_default"),
    n_max::Union{Nothing, Int}=nothing,
    p_num::Union{Nothing, Int}=nothing,
    pz_max::Union{Nothing, Real}=nothing,
    cutoff_N::Union{Nothing, Int}=nothing,
    route::Union{Nothing, Symbol}=nothing,
    zeta_num::Union{Nothing, Int}=nothing,
    n_max_policy::Union{Nothing, Symbol}=nothing,
    thermal_tail_factor::Union{Nothing, Real}=nothing,
    n_max_floor::Union{Nothing, Int}=nothing,
    n_max_cap::Union{Nothing, Int}=nothing,
    imc=nothing,
    kwargs...,
)
    isempty(kwargs) || throw(ArgumentError("unsupported PNJLMagneticModel keyword(s): $(collect(keys(kwargs)))"))
    base = PNJLModel(; profile=profile, physics_profile=physics_profile)
    thermo = _magnetic_thermodynamics_module()
    config_defaults = thermo.default_magnetic_config(
        eB_fm2=float(eB_fm2),
        profile=magnetic_profile,
        params=base.params,
    )
    imc_params = imc === nothing ? config_defaults.imc : imc
    conf = thermo.MagneticConfig(
        eB_fm2=float(eB_fm2),
        n_max=n_max === nothing ? config_defaults.n_max : n_max,
        p_num=p_num === nothing ? config_defaults.p_num : p_num,
        pz_max=pz_max === nothing ? config_defaults.pz_max : float(pz_max),
        cutoff_N=cutoff_N === nothing ? config_defaults.cutoff_N : cutoff_N,
        imc=imc_params,
        route=route === nothing ? config_defaults.route : route,
        zeta_num=zeta_num === nothing ? config_defaults.zeta_num : zeta_num,
        n_max_policy=n_max_policy === nothing ? config_defaults.n_max_policy : n_max_policy,
        thermal_tail_factor=thermal_tail_factor === nothing ? config_defaults.thermal_tail_factor : thermal_tail_factor,
        n_max_floor=n_max_floor === nothing ? config_defaults.n_max_floor : n_max_floor,
        n_max_cap=n_max_cap === nothing ? config_defaults.n_max_cap : n_max_cap,
        params=base.params,
    )
    return PNJLMagneticModel(base, conf)
end

@inline function calculate_mass_vec(model::PNJLMagneticModel, φ; kwargs...)
    thermo = _magnetic_thermodynamics_module()
    # Preserve Dual/other Real element types for generic derivative callers.
    φ3 = SVector(φ[1], φ[2], φ[3])
    G_B = thermo.coupling_GB(model.magnetic.eB_fm2; G0=model.base.params.G_fm2, imc=model.magnetic.imc)
    return thermo._calculate_mass_vec_with_GB(φ3, G_B; params=model.base.params)
end

@inline function calculate_chiral(model::PNJLMagneticModel, φ; kwargs...)
    thermo = _magnetic_thermodynamics_module()
    φ3 = SVector(φ[1], φ[2], φ[3])
    G_B = thermo.coupling_GB(model.magnetic.eB_fm2; G0=model.base.params.G_fm2, imc=model.magnetic.imc)
    return thermo._chiral_with_GB(φ3, G_B; params=model.base.params)
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

@inline function number_densities(model::PNJLMagneticModel, x_state, T, mu_vec; thermal_nodes=nothing, kwargs...)
    _ = thermal_nodes
    return _magnetic_thermodynamics_module().calculate_magnetic_number_densities(
        x_state,
        normalize_mu_vec(mu_vec),
        T,
        model.magnetic;
        kwargs...,
    )
end

function model_rho(
    model::PNJLMagneticModel,
    x_state,
    mu_vec,
    T_fm;
    p_num::Union{Nothing, Int}=nothing,
    t_num::Int=8,
    xi::Real=0.0,
    thermal_nodes=nothing,
    kwargs...,
)
    _ = thermal_nodes
    return _magnetic_thermodynamics_module().calculate_magnetic_rho(
        x_state,
        normalize_mu_vec(mu_vec),
        T_fm,
        model.magnetic;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        kwargs...,
    )
end

function model_thermo(
    model::PNJLMagneticModel,
    x_state,
    mu_vec,
    T_fm;
    p_num::Union{Nothing, Int}=nothing,
    t_num::Int=8,
    xi::Real=0.0,
    thermal_nodes=nothing,
    entropy_step::Real=1e-4,
    kwargs...,
)
    _ = thermal_nodes
    T_value = Float64(T_fm)
    isfinite(T_value) && T_value > 0.0 || throw(ArgumentError(
        "magnetic thermodynamics requires finite T_fm > 0, got $(T_fm)",
    ))
    h = Float64(entropy_step)
    isfinite(h) && h > 0.0 || throw(ArgumentError("entropy_step must be positive"))
    μ = normalize_mu_vec(mu_vec)
    pressure = -omega(model, x_state, T_fm, μ; p_num=p_num, t_num=t_num, xi=xi, kwargs...)
    rho_vec = model_rho(model, x_state, μ, T_fm; p_num=p_num, t_num=t_num, xi=xi, kwargs...)
    Tp = T_value + h
    Tm = max(T_value - h, h / 10)
    entropy = (
        -omega(model, x_state, Tp, μ; p_num=p_num, t_num=t_num, xi=xi, kwargs...) -
        (-omega(model, x_state, Tm, μ; p_num=p_num, t_num=t_num, xi=xi, kwargs...))
    ) / (Tp - Tm)
    rho_norm = sum(rho_vec) / (3.0 * Main.Constants_PNJL.ρ0_inv_fm3)
    energy = -pressure + sum(μ .* rho_vec) + T_value * entropy
    return pressure, rho_norm, entropy, energy
end

function omega_components(model::PNJLMagneticModel, x_state, T, mu_vec; kwargs...)
    st = x_state isa MeanFieldState ? state_vector(x_state) : SVector{5, Float64}(Tuple(x_state))
    μ = normalize_mu_vec(mu_vec)
    comp = _magnetic_thermodynamics_module().calculate_magnetic_omega_components(st, μ, T, model.magnetic; kwargs...)
    return (
        chi=comp.chi,
        poly=comp.poly,
        vac=comp.vac,
        therm=comp.therm,
        masses=comp.masses,
        omega=comp.omega,
        n_max=comp.n_max,
        G_B=comp.G_B,
        route=model.magnetic.route,
        zeta_num=model.magnetic.zeta_num,
    )
end

@inline function omega(model::PNJLMagneticModel, x_state, T, mu_vec; kwargs...)
    return omega_components(model, x_state, T, mu_vec; kwargs...).omega
end
