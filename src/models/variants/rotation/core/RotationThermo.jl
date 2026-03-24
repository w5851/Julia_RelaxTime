module RotationThermo

using SpecialFunctions
using StaticArrays

const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "config", "ConfigLoader.jl"))
if !isdefined(@__MODULE__, :ConfigLoader)
    include(_CONFIG_LOADER_PATH)
end
using .ConfigLoader: load_config

export RotationParams
export effective_energy_rot
export pressure_density_entropy_energy
export omega_components
export pressure_derivative_omega
export gap_residuals
export polyakov_potential

const _DEFAULT_PROFILE = "default"
const _ROT_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "..", "config", "models", "rotation"))
const _DEFAULT_ROT_CONFIG = Dict{String, Any}(
    "model" => Dict{String, Any}(
        "m0_MeV" => 5.0,
        "Lambda_MeV" => 650.0,
        "G_MeV_inv2" => 4.93,
        "a0" => 6.75,
        "a1" => -1.95,
        "a2" => 2.625,
        "a3" => -7.44,
        "b3" => 0.75,
        "b4" => 7.5,
        "T0_MeV" => 270.0,
        "r_fm" => 0.1,
        "n_cut" => 5,
        "hbarc_MeV_fm" => 197.3269804,
    ),
)

@inline _mev_to_fminv(x_MeV::Real, hbarc_MeV_fm::Real) = float(x_MeV) / float(hbarc_MeV_fm)

@inline function _load_rotation_config(; profile::String=get(ENV, "ROTATION_PARAM_PROFILE", _DEFAULT_PROFILE))
    return load_config(_ROT_CONFIG_DIR, _DEFAULT_ROT_CONFIG; profile=profile).config
end

Base.@kwdef struct RotationParams
    hbarc_MeV_fm::Float64 = 197.3269804
    m0_inv_fm::Float64 = 5.0 / 197.3269804
    Λ_inv_fm::Float64 = 650.0 / 197.3269804
    G_fm2::Float64 = 4.93 * 197.3269804^2

    a0::Float64 = 6.75
    a1::Float64 = -1.95
    a2::Float64 = 2.625
    a3::Float64 = -7.44
    b3::Float64 = 0.75
    b4::Float64 = 7.5
    T0_inv_fm::Float64 = 270.0 / 197.3269804

    r_fm::Float64 = 0.1
    n_cut::Int = 5
end

function RotationParams(; profile::Union{Nothing, String}=nothing, kwargs...)
    profile_eff = something(profile, get(ENV, "ROTATION_PARAM_PROFILE", _DEFAULT_PROFILE))
    model = get(_load_rotation_config(profile=profile_eff), "model", Dict{String, Any}())
    hbarc = Float64(get(model, "hbarc_MeV_fm", 197.3269804))

    base = (
        hbarc_MeV_fm=hbarc,
        m0_inv_fm=_mev_to_fminv(Float64(get(model, "m0_MeV", 5.0)), hbarc),
        Λ_inv_fm=_mev_to_fminv(Float64(get(model, "Lambda_MeV", 650.0)), hbarc),
        G_fm2=Float64(get(model, "G_MeV_inv2", 4.93)) * hbarc^2,
        a0=Float64(get(model, "a0", 6.75)),
        a1=Float64(get(model, "a1", -1.95)),
        a2=Float64(get(model, "a2", 2.625)),
        a3=Float64(get(model, "a3", -7.44)),
        b3=Float64(get(model, "b3", 0.75)),
        b4=Float64(get(model, "b4", 7.5)),
        T0_inv_fm=_mev_to_fminv(Float64(get(model, "T0_MeV", 270.0)), hbarc),
        r_fm=Float64(get(model, "r_fm", 0.1)),
        n_cut=Int(get(model, "n_cut", 5)),
    )
    merged = merge(base, kwargs)
    return RotationParams(
        Float64(merged[:hbarc_MeV_fm]),
        Float64(merged[:m0_inv_fm]),
        Float64(merged[:Λ_inv_fm]),
        Float64(merged[:G_fm2]),
        Float64(merged[:a0]),
        Float64(merged[:a1]),
        Float64(merged[:a2]),
        Float64(merged[:a3]),
        Float64(merged[:b3]),
        Float64(merged[:b4]),
        Float64(merged[:T0_inv_fm]),
        Float64(merged[:r_fm]),
        Int(merged[:n_cut]),
    )
end

@inline function effective_energy_rot(p_t::Real, p_z::Real, mass::Real, mode_n::Int, omega::Real)
    return sqrt(p_t^2 + p_z^2 + mass^2) - (mode_n + 0.5) * omega
end

@inline function polyakov_potential(Φ::Real, Φbar::Real, T::Real, p::RotationParams)
    T <= eps(Float64) && return 0.0
    x = p.T0_inv_fm / T
    b2 = p.a0 + p.a1 * x + p.a2 * x^2 + p.a3 * x^3
    return T^4 * (-0.5 * b2 * Φ * Φbar - (p.b3 / 6) * (Φ^3 + Φbar^3) + 0.25 * p.b4 * (Φ * Φbar)^2)
end

@inline function _qplus(E::Real, μ::Real, T::Real, Φ::Real, Φbar::Real)
    ex = exp(-(E - μ) / T)
    return log1p(3 * Φ * ex + 3 * Φbar * ex^2 + ex^3)
end

@inline function _qminus(E::Real, μ::Real, T::Real, Φ::Real, Φbar::Real)
    ex = exp(-(E + μ) / T)
    return log1p(3 * Φbar * ex + 3 * Φ * ex^2 + ex^3)
end

function _thermal_sum(mass::Real, T::Real, μ::Real, omega::Real, Φ::Real, Φbar::Real, p::RotationParams; p_num::Int=48, z_num::Int=48)
    T <= eps(Float64) && return 0.0

    pt_max = p.Λ_inv_fm
    dpt = pt_max / p_num
    total = 0.0

    @inbounds for n in -p.n_cut:p.n_cut
        for i in 0:p_num
            pt = i * dpt
            wpt = (i == 0 || i == p_num) ? 0.5 : 1.0
            zmax = sqrt(max(p.Λ_inv_fm^2 - pt^2, 0.0))
            if zmax <= 0
                continue
            end
            dz = 2zmax / z_num
            wn = besselj(n, pt * p.r_fm)^2 + besselj(n + 1, pt * p.r_fm)^2

            for j in 0:z_num
                pz = -zmax + j * dz
                wz = (j == 0 || j == z_num) ? 0.5 : 1.0
                e = effective_energy_rot(pt, pz, mass, n, omega)
                total += wpt * wz * pt * wn * (_qplus(e, μ, T, Φ, Φbar) + _qminus(e, μ, T, Φ, Φbar))
            end
        end
    end
    return (T / (4pi^2)) * dpt * (2 * p.Λ_inv_fm / z_num) * total
end

@inline function omega_components(phi::Real, Φ::Real, Φbar::Real, T::Real, μ::Real, omega::Real, p::RotationParams; p_num::Int=48, z_num::Int=48)
    mass = p.m0_inv_fm - 2 * p.G_fm2 * phi
    chi = p.G_fm2 * phi^2
    poly = polyakov_potential(Φ, Φbar, T, p)
    therm = -_thermal_sum(mass, T, μ, omega, Φ, Φbar, p; p_num=p_num, z_num=z_num)
    vac = 0.0
    omega_total = chi + poly + vac + therm
    Tm = typeof(mass)
    masses = SVector{3, Tm}(mass, mass, mass)
    return (chi=chi, vac=vac, therm=therm, poly=poly, masses=masses, omega=omega_total)
end

@inline function omega_components(phi::Real, T::Real, μ::Real, omega::Real, p::RotationParams; p_num::Int=48, z_num::Int=48)
    return omega_components(phi, 1.0, 1.0, T, μ, omega, p; p_num=p_num, z_num=z_num)
end

@inline function pressure_density_entropy_energy(phi::Real, T::Real, μ::Real, omega::Real, p::RotationParams; Φ::Real=1.0, Φbar::Real=1.0, p_num::Int=48, z_num::Int=48)
    comp = omega_components(phi, Φ, Φbar, T, μ, omega, p; p_num=p_num, z_num=z_num)
    pressure = -comp.omega

    hμ = max(1e-6, 1e-3 * max(abs(μ), 1.0))
    p_plus = -omega_components(phi, Φ, Φbar, T, μ + hμ, omega, p; p_num=p_num, z_num=z_num).omega
    p_minus = -omega_components(phi, Φ, Φbar, T, μ - hμ, omega, p; p_num=p_num, z_num=z_num).omega
    rho = (p_plus - p_minus) / (2hμ)

    hT = max(1e-6, 1e-3 * max(T, 1.0))
    t_plus = -omega_components(phi, Φ, Φbar, T + hT, μ, omega, p; p_num=p_num, z_num=z_num).omega
    t_minus = -omega_components(phi, Φ, Φbar, max(T - hT, 1e-8), μ, omega, p; p_num=p_num, z_num=z_num).omega
    entropy = (t_plus - t_minus) / (2hT)

    energy = -pressure + T * entropy + μ * rho
    return (pressure=pressure, rho=rho, entropy=entropy, energy=energy, mass=comp.masses[1])
end

@inline function gap_residuals(phi::Real, Φ::Real, Φbar::Real, T::Real, μ::Real, omega::Real, p::RotationParams; p_num::Int=48, z_num::Int=48)
    h = 1e-4
    dphi = (
        omega_components(phi + h, Φ, Φbar, T, μ, omega, p; p_num=p_num, z_num=z_num).omega -
        omega_components(phi - h, Φ, Φbar, T, μ, omega, p; p_num=p_num, z_num=z_num).omega
    ) / (2h)
    dPhi = (
        omega_components(phi, Φ + h, Φbar, T, μ, omega, p; p_num=p_num, z_num=z_num).omega -
        omega_components(phi, Φ - h, Φbar, T, μ, omega, p; p_num=p_num, z_num=z_num).omega
    ) / (2h)
    dPhiBar = (
        omega_components(phi, Φ, Φbar + h, T, μ, omega, p; p_num=p_num, z_num=z_num).omega -
        omega_components(phi, Φ, Φbar - h, T, μ, omega, p; p_num=p_num, z_num=z_num).omega
    ) / (2h)
    Tr = promote_type(typeof(dphi), typeof(dPhi), typeof(dPhiBar))
    return SVector{3, Tr}(dphi, dPhi, dPhiBar)
end

@inline function pressure_derivative_omega(phi::Real, Φ::Real, Φbar::Real, T::Real, μ::Real, omega::Real, p::RotationParams)
    h = 1e-4
    p_plus = -omega_components(phi, Φ, Φbar, T, μ, omega + h, p).omega
    p_minus = -omega_components(phi, Φ, Φbar, T, μ, omega - h, p).omega
    return (p_plus - p_minus) / (2h)
end

@inline function pressure_derivative_omega(phi::Real, T::Real, μ::Real, omega::Real, p::RotationParams)
    return pressure_derivative_omega(phi, 1.0, 1.0, T, μ, omega, p)
end

end # module
