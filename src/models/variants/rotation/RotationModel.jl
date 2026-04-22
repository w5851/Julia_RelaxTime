"""RotationModel

Rotation-PNJL 适配壳（最小可运行版）。
"""

using NLsolve
using StaticArrays

const _ROTATION_THERMO_PATH = normpath(joinpath(@__DIR__, "core", "RotationThermo.jl"))
if !isdefined(@__MODULE__, :RotationThermo)
    include(_ROTATION_THERMO_PATH)
end

using .RotationThermo: RotationParams
using .RotationThermo: pressure_density_entropy_energy as rotation_pressure_density_entropy_energy
using .RotationThermo: omega_components as rotation_omega_components
using .RotationThermo: gap_residuals as rotation_gap_residuals
using .RotationThermo: polyakov_potential as rotation_polyakov_potential

export RotationModel

struct RotationModel <: AbstractQCDModel
    params::RotationParams
end

@inline model_capabilities(::RotationModel) = ModelCapabilities(
    supports_solve_gap=true,
    supports_model_thermo=true,
    supports_number_densities=true,
)

RotationModel(; kwargs...) = RotationModel(RotationParams(; kwargs...))

@inline gap_state_dim(::RotationModel) = 3

function solve_gap(model::RotationModel, T_fm, mu_vec; omega::Real=0.0, initial_guess=nothing, kwargs...)
    p_num = Int(get(kwargs, :p_num, 48))
    z_num = Int(get(kwargs, :z_num, 48))
    mu = mu_vec isa Real ? mu_vec : (mu_vec[1] + mu_vec[2] + mu_vec[3]) / 3
    x0 = initial_guess === nothing ? [0.1, 1.0, 1.0] : collect(float.(initial_guess))

    function residual!(F, x)
        r = rotation_gap_residuals(x[1], x[2], x[3], T_fm, mu, omega, model.params; p_num=p_num, z_num=z_num)
        F[1] = r[1]
        F[2] = r[2]
        F[3] = r[3]
        return nothing
    end

    sol = nlsolve(residual!, x0; autodiff=:forward, method=:newton, xtol=1e-10, ftol=1e-10)
    if !(sol.f_converged && isfinite(sol.residual_norm))
        error("rotation solve_gap failed: f_converged=$(sol.f_converged), residual_norm=$(sol.residual_norm)")
    end

    x = sol.zero
    phi = SVector{3, Float64}(float(x[1]), 0.0, 0.0)
    return MeanFieldState(phi; Phi=float(x[2]), PhiBar=float(x[3]))
end

@inline function _mu_scalar(mu_vec)
    return mu_vec isa Real ? mu_vec : (mu_vec[1] + mu_vec[2] + mu_vec[3]) / 3
end

function calculate_mass_vec(model::RotationModel, φ; kwargs...)
    _ = kwargs
    m = model.params.m0_inv_fm - 2 * model.params.G_fm2 * φ[1]
    return SVector{3, typeof(m)}(m, m, m)
end

@inline calculate_chiral(model::RotationModel, φ; kwargs...) = model.params.G_fm2 * (φ[1]^2)

@inline function polyakov_potential(model::RotationModel, Φ, Φbar, T; kwargs...)
    _ = kwargs
    return rotation_polyakov_potential(Φ, Φbar, T, model.params)
end

@inline function vacuum_contribution(::RotationModel, masses; kwargs...)
    _ = (masses, kwargs)
    return 0.0
end

@inline function _phi_from_masses(model::RotationModel, masses)
    return (model.params.m0_inv_fm - masses[1]) / (2 * model.params.G_fm2)
end

function thermal_contribution(model::RotationModel, masses, Φ, Φbar, mu_vec, T; omega::Real=0.0, kwargs...)
    p_num = Int(get(kwargs, :p_num, 48))
    z_num = Int(get(kwargs, :z_num, 48))
    mu = _mu_scalar(mu_vec)
    phi = _phi_from_masses(model, masses)
    return rotation_omega_components(phi, Φ, Φbar, T, mu, omega, model.params; p_num=p_num, z_num=z_num).therm
end

function number_densities(model::RotationModel, x_state, T, mu_vec; omega::Real=0.0, kwargs...)
    p_num = Int(get(kwargs, :p_num, 48))
    z_num = Int(get(kwargs, :z_num, 48))
    st = x_state isa MeanFieldState ? x_state : MeanFieldState(x_state)
    mu = _mu_scalar(mu_vec)
    thermo = rotation_pressure_density_entropy_energy(st.phi[1], T, mu, omega, model.params; Φ=st.Phi, Φbar=st.PhiBar, p_num=p_num, z_num=z_num)
    Tm = typeof(thermo.rho)
    q = SVector{3, Tm}(thermo.rho, thermo.rho, thermo.rho)
    aq = SVector{3, Tm}(zero(Tm), zero(Tm), zero(Tm))
    return (quark=q, antiquark=aq)
end
