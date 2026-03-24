module GasLiquidThermodynamics

using StaticArrays
using ..GasLiquidEquationSet: GasLiquidCoreParams, GasLiquidState, effective_masses, density_bundle

export effective_masses
export pressure_density_entropy_energy
export omega_components

@inline function _scalar_potential(st::GasLiquidState, p::GasLiquidCoreParams)
    x = p.g_sigma * st.sigma
    return (1 / 3) * p.b * p.m_nucleon_inv_fm * x^3 + (1 / 4) * p.c * x^4
end

@inline function _vector_fields(st::GasLiquidState, p::GasLiquidCoreParams)
    omega0 = (st.mu_p + st.mu_n) / (2 * p.g_omega + eps(Float64))
    rho03 = (st.mu_n - st.mu_p) / (2 * p.g_rho + eps(Float64))
    return omega0, rho03
end

@inline function pressure_density_entropy_energy(st::GasLiquidState, T::Real, p::GasLiquidCoreParams; p_num::Int=96)
    dens = density_bundle(st, T, p; p_num=p_num)
    Tv = float(T)

    omega0, rho03 = _vector_fields(st, p)
    us = _scalar_potential(st, p)

    omega_val = us - 0.5 * p.m_omega_inv_fm^2 * omega0^2 - 0.5 * p.m_rho_inv_fm^2 * rho03^2 - 0.5 * p.m_delta_inv_fm^2 * st.delta^2
    if Tv > eps(Float64)
        omega_val -= 2.0 * Tv * dens.logsum
    end

    pressure = -omega_val
    rhoB = dens.rho_B

    # finite-difference entropy at fixed state/mu (for thermo_kernel compatibility)
    h = max(1e-5, 1e-3 * max(Tv, 1e-3))
    dens_p = density_bundle(st, Tv + h, p; p_num=p_num)
    dens_m = density_bundle(st, max(Tv - h, 1e-8), p; p_num=p_num)

    op = us - 0.5 * p.m_omega_inv_fm^2 * omega0^2 - 0.5 * p.m_rho_inv_fm^2 * rho03^2 - 0.5 * p.m_delta_inv_fm^2 * st.delta^2 - 2.0 * (Tv + h) * dens_p.logsum
    om = us - 0.5 * p.m_omega_inv_fm^2 * omega0^2 - 0.5 * p.m_rho_inv_fm^2 * rho03^2 - 0.5 * p.m_delta_inv_fm^2 * st.delta^2 - 2.0 * max(Tv - h, 1e-8) * dens_m.logsum
    entropy = -(op - om) / (2h)

    energy = us + 0.5 * p.m_omega_inv_fm^2 * omega0^2 + 0.5 * p.m_rho_inv_fm^2 * rho03^2 + 0.5 * p.m_delta_inv_fm^2 * st.delta^2 + dens.kinetic_energy

    Tm = promote_type(typeof(dens.masses.mp), typeof(dens.masses.mn), typeof(p.m_nucleon_inv_fm))
    masses = SVector{3, Tm}(dens.masses.mp, dens.masses.mn, p.m_nucleon_inv_fm)
    return (pressure=pressure, rho=rhoB, entropy=entropy, energy=energy, masses=masses)
end

@inline function omega_components(st::GasLiquidState, T::Real, p::GasLiquidCoreParams; p_num::Int=96)
    thermo = pressure_density_entropy_energy(st, T, p; p_num=p_num)
    omega0, rho03 = _vector_fields(st, p)

    chi = _scalar_potential(st, p) + 0.5 * p.m_delta_inv_fm^2 * st.delta^2
    vac = -0.5 * p.m_omega_inv_fm^2 * omega0^2 - 0.5 * p.m_rho_inv_fm^2 * rho03^2
    therm = -thermo.pressure - chi - vac
    poly = 0.0
    omega = chi + vac + therm + poly
    return (chi=chi, vac=vac, therm=therm, poly=poly, masses=thermo.masses, omega=omega)
end

end # module
