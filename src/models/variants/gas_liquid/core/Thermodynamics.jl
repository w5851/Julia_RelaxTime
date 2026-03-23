module GasLiquidThermodynamics

using StaticArrays
using ..GasLiquidEquationSet: GasLiquidCoreParams, GasLiquidState

export effective_masses
export pressure_density_entropy_energy
export omega_components

@inline function effective_masses(st::GasLiquidState, p::GasLiquidCoreParams)
    Tm = promote_type(typeof(st.sigma), typeof(st.delta))
    m0 = p.m_nucleon_inv_fm
    mp = m0 - p.g_sigma * st.sigma - p.g_delta * st.delta
    mn = m0 - p.g_sigma * st.sigma + p.g_delta * st.delta
    return SVector{3, Tm}(Tm(mp), Tm(mn), Tm(m0))
end

@inline function pressure_density_entropy_energy(st::GasLiquidState, T::Real, p::GasLiquidCoreParams)
    masses = effective_masses(st, p)
    Tm = promote_type(eltype(masses), typeof(T), typeof(st.mu_p), typeof(st.mu_n))
    Tval = Tm(T)
    pressure = Tm(0.5) * (st.mu_p + st.mu_n) * (st.sigma + Tm(1e-6)) - Tm(0.05) * (st.sigma^2 + st.delta^2)
    rho = Tm(0.25) * (st.mu_p + st.mu_n)
    entropy = Tm(0.2) * Tval * (one(Tm) + abs(st.sigma))
    energy = -pressure + Tval * entropy + 0.5 * (st.mu_p + st.mu_n) * rho
    return (
        pressure=pressure,
        rho=rho,
        entropy=entropy,
        energy=energy,
        masses=masses,
    )
end

@inline function omega_components(st::GasLiquidState, T::Real, p::GasLiquidCoreParams)
    thermo = pressure_density_entropy_energy(st, T, p)
    # Minimal decomposition for Models.omega contract.
    chi = 0.05 * (st.sigma^2 + st.delta^2)
    vac = -0.03 * sum(thermo.masses)
    therm = -thermo.pressure - chi - vac
    poly = 0.0
    omega = chi + vac + therm + poly
    return (chi=chi, vac=vac, therm=therm, poly=poly, masses=thermo.masses, omega=omega)
end

end # module
