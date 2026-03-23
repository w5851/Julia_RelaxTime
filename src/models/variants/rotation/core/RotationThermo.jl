module RotationThermo

using StaticArrays

export RotationParams
export effective_energy_rot
export pressure_density_entropy_energy
export omega_components
export pressure_derivative_omega

Base.@kwdef struct RotationParams
    m0_inv_fm::Float64 = 300.0 / 197.3269804
    g_chiral::Float64 = 1.8
end

@inline function effective_energy_rot(p::Real, mass::Real, mode_n::Int, omega::Real)
    return sqrt(p^2 + mass^2) - (mode_n + 0.5) * omega
end

@inline function pressure_density_entropy_energy(phi::Real, T::Real, mu::Real, omega::Real, params::RotationParams)
    mass = params.m0_inv_fm - params.g_chiral * phi
    e0 = effective_energy_rot(0.2, mass, 0, omega)
    pressure = 0.25 * mu * (phi + 1e-6) - 0.05 * phi^2 + 0.02 * omega * e0
    rho = 0.30 * mu * (1 + 0.1 * omega)
    entropy = 0.15 * T * (1 + abs(phi))
    energy = -pressure + T * entropy + mu * rho
    return (pressure=pressure, rho=rho, entropy=entropy, energy=energy, mass=mass)
end

@inline function omega_components(phi::Real, T::Real, mu::Real, omega::Real, params::RotationParams)
    thermo = pressure_density_entropy_energy(phi, T, mu, omega, params)
    chi = 0.05 * phi^2
    vac = -0.04 * thermo.mass
    therm = -thermo.pressure - chi - vac
    poly = 0.0
    omega_total = chi + vac + therm + poly
    masses = SVector{3, typeof(thermo.mass)}(thermo.mass, thermo.mass, params.m0_inv_fm)
    return (chi=chi, vac=vac, therm=therm, poly=poly, masses=masses, omega=omega_total)
end

@inline function pressure_derivative_omega(phi::Real, T::Real, mu::Real, omega::Real, params::RotationParams)
    h = 1e-4
    p_plus = pressure_density_entropy_energy(phi, T, mu, omega + h, params).pressure
    p_minus = pressure_density_entropy_energy(phi, T, mu, omega - h, params).pressure
    return (p_plus - p_minus) / (2h)
end

end # module
