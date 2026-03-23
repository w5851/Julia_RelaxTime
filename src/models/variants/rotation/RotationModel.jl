"""RotationModel

Rotation-PNJL 适配壳（最小可运行版）。
"""

using NLsolve
using StaticArrays

const _ROTATION_THERMO_PATH = normpath(joinpath(@__DIR__, "core", "RotationThermo.jl"))
if !isdefined(@__MODULE__, :RotationThermo)
    include(_ROTATION_THERMO_PATH)
end

using .RotationThermo: RotationParams, pressure_density_entropy_energy, omega_components as rotation_omega_components

export RotationModel

struct RotationModel <: AbstractQCDModel
    params::RotationParams
end

RotationModel(; kwargs...) = RotationModel(RotationParams(; kwargs...))

@inline gap_state_dim(::RotationModel) = 3

function solve_gap(model::RotationModel, T_fm, mu_vec; omega::Real=0.0, initial_guess=nothing, kwargs...)
    _ = (kwargs,)
    mu = mu_vec isa Real ? mu_vec : (mu_vec[1] + mu_vec[2] + mu_vec[3]) / 3
    x0 = initial_guess === nothing ? [0.1, 1.0, 1.0] : collect(float.(initial_guess))

    function residual!(F, x)
        F[1] = x[1] - tanh((mu - T_fm) / (T_fm + 0.2))
        F[2] = x[2] - 1.0
        F[3] = x[3] - 1.0
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
    _ = (kwargs,)
    m = model.params.m0_inv_fm - model.params.g_chiral * φ[1]
    return SVector{3, typeof(m)}(m, m, model.params.m0_inv_fm)
end

@inline calculate_chiral(::RotationModel, φ; kwargs...) = 0.05 * (φ[1]^2)

@inline function polyakov_potential(::RotationModel, Φ, Φbar, T; kwargs...)
    _ = (Φ, Φbar, T, kwargs)
    return 0.0
end

@inline vacuum_contribution(::RotationModel, masses; kwargs...) = -0.04 * masses[1]

function thermal_contribution(model::RotationModel, masses, Φ, Φbar, mu_vec, T; omega::Real=0.0, kwargs...)
    _ = (masses, Φ, Φbar, kwargs)
    mu = _mu_scalar(mu_vec)
    thermo = pressure_density_entropy_energy(0.0, T, mu, omega, model.params)
    chi = 0.0
    vac = -0.04 * model.params.m0_inv_fm
    return -thermo.pressure - chi - vac
end

function number_densities(model::RotationModel, x_state, T, mu_vec; omega::Real=0.0, kwargs...)
    _ = (kwargs,)
    st = x_state isa MeanFieldState ? x_state : MeanFieldState(x_state)
    mu = _mu_scalar(mu_vec)
    thermo = pressure_density_entropy_energy(st.phi[1], T, mu, omega, model.params)
    Tm = typeof(thermo.rho)
    q = SVector{3, Tm}(thermo.rho, thermo.rho, thermo.rho)
    aq = SVector{3, Tm}(zero(Tm), zero(Tm), zero(Tm))
    return (quark=q, antiquark=aq)
end

function omega_components(model::RotationModel, x_state, T, mu_vec; omega::Real=0.0, kwargs...)
    _ = (kwargs,)
    st = x_state isa MeanFieldState ? x_state : MeanFieldState(x_state)
    mu = _mu_scalar(mu_vec)
    return rotation_omega_components(st.phi[1], T, mu, omega, model.params)
end
