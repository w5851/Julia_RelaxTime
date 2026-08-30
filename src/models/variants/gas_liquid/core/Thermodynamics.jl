"""Thermodynamic observables derived from the RMF grand potential."""
module GasLiquidThermodynamics

using StaticArrays
using ..GasLiquidEquationSet: GasLiquidCoreParams, GasLiquidState, density_bundle
using ..GasLiquidEquationSet: effective_masses, field_contributions, physical_mu_pair

export effective_masses
export pressure_density_entropy_energy, omega_components
export scalar_potential, scalar_field_contribution, omega_identity_residual
export thermodynamic_consistency_report

@inline function scalar_potential(st::GasLiquidState, p::GasLiquidCoreParams)
    # b/c follow the thesis convention U_sigma = b*m_N*S^3/3 + c*S^4/4;
    # the powers of g_sigma are therefore not absorbed into the coefficients.
    x = st.S
    return (p.b * p.m_nucleon_inv_fm * x^3) / 3.0 + (p.c * x^4) / 4.0
end

@inline function _delta_potential(st::GasLiquidState, p::GasLiquidCoreParams)
    p.f_delta == 0.0 ? 0.0 : -st.D^2 / (2.0 * p.f_delta)
end

@inline function _scalar_field_contribution(st::GasLiquidState, p::GasLiquidCoreParams)
    scalar_potential(st, p) - (p.f_sigma == 0.0 ? 0.0 : st.S^2 / (2.0 * p.f_sigma)) + _delta_potential(st, p)
end

@inline scalar_field_contribution(st::GasLiquidState, p::GasLiquidCoreParams) = _scalar_field_contribution(st, p)

@inline function _vector_field_contribution(dens, p::GasLiquidCoreParams)
    fields = field_contributions(dens, p)
    omega = p.f_omega == 0.0 ? 0.0 : -fields.W^2 / (2.0 * p.f_omega)
    rho = p.f_rho == 0.0 ? 0.0 : -fields.R^2 / (2.0 * p.f_rho)
    return (omega=omega, rho=rho, W=fields.W, R=fields.R)
end

@inline function _omega_fixed_fields(st::GasLiquidState, T::Real, p::GasLiquidCoreParams, W::Real, R::Real;
    p_num::Int=96, p_max_inv_fm::Real=10.0)
    dens = density_bundle(st, T, p; p_num=p_num, p_max_inv_fm=p_max_inv_fm)
    scalar = _scalar_field_contribution(st, p)
    vector = p.f_omega == 0.0 ? 0.0 : -(Float64(W)^2 / (2.0 * p.f_omega))
    vector += p.f_rho == 0.0 ? 0.0 : -(Float64(R)^2 / (2.0 * p.f_rho))
    return scalar + vector - 2.0 * float(T) * dens.thermal_log_integral
end

"""Evaluate the common Omega/P/rho/s/epsilon contract.

The entropy contains both particle and antiparticle terms. `energy` is computed
from `epsilon + P = T*s + mu_p*rho_p + mu_n*rho_n`, so the same convention is
used for fixed-mu and fixed-rho workflows.
"""
function pressure_density_entropy_energy(st::GasLiquidState, T::Real, p::GasLiquidCoreParams;
    p_num::Int=96, p_max_inv_fm::Real=10.0,
    mu_p_target::Union{Nothing, Real}=nothing,
    mu_n_target::Union{Nothing, Real}=nothing)
    dens = density_bundle(st, T, p; p_num=p_num, p_max_inv_fm=p_max_inv_fm)
    scalar = _scalar_field_contribution(st, p)
    vector = _vector_field_contribution(dens, p)
    Tval = float(T)
    thermal = -2.0 * Tval * dens.thermal_log_integral
    omega = scalar + vector.omega + vector.rho + thermal
    pressure = -omega

    mu = if mu_p_target === nothing && mu_n_target === nothing
        physical_mu_pair(st, dens, p)
    elseif mu_p_target !== nothing && mu_n_target !== nothing
        (mu_p=Float64(mu_p_target), mu_n=Float64(mu_n_target))
    else
        throw(ArgumentError("mu_p_target and mu_n_target must be supplied together"))
    end
    entropy = dens.entropy
    energy = Tval * entropy + mu.mu_p * dens.rho_p + mu.mu_n * dens.rho_n - pressure

    Tm = promote_type(typeof(dens.masses.mp), typeof(dens.masses.mn), typeof(p.m_nucleon_inv_fm))
    masses = SVector{3, Tm}(dens.masses.mp, dens.masses.mn, p.m_nucleon_inv_fm)
    return (
        omega=omega,
        pressure=pressure,
        rho=dens.rho_B,
        rho_p=dens.rho_p,
        rho_n=dens.rho_n,
        rho_B=dens.rho_B,
        rho_3=dens.rho_3,
        rho_s_p=dens.rho_s_p,
        rho_s_n=dens.rho_s_n,
        rho_s3=dens.rho_s3,
        entropy=entropy,
        energy=energy,
        kinetic_energy=dens.kinetic_energy,
        masses=masses,
        mu_p=mu.mu_p,
        mu_n=mu.mu_n,
        mu_p_star=st.mu_p,
        mu_n_star=st.mu_n,
        S=st.S,
        D=st.D,
        W=vector.W,
        R=vector.R,
        chi=scalar,
        vac=vector.omega + vector.rho,
        therm=thermal,
        quadrature=(p_num=p_num, p_max_inv_fm=Float64(p_max_inv_fm)),
    )
end

"""Return the decomposed Omega used by the Models compatibility layer."""
@inline function omega_components(st::GasLiquidState, T::Real, p::GasLiquidCoreParams;
    p_num::Int=96, p_max_inv_fm::Real=10.0)
    thermo = pressure_density_entropy_energy(st, T, p; p_num=p_num, p_max_inv_fm=p_max_inv_fm)
    return (
        chi=thermo.chi,
        poly=0.0,
        # The generic Models adapter may call vacuum_contribution without T.
        # Keep that compatibility component T-independent and retain the
        # density-dependent vector contribution in the thermal aggregate.
        vac=0.0,
        therm=thermo.vac + thermo.therm,
        vector_vac=thermo.vac,
        masses=thermo.masses,
        omega=thermo.omega,
        pressure=thermo.pressure,
        rho=thermo.rho,
        rho_p=thermo.rho_p,
        rho_n=thermo.rho_n,
        rho_3=thermo.rho_3,
        rho_s3=thermo.rho_s3,
        entropy=thermo.entropy,
        energy=thermo.energy,
        S=thermo.S,
        D=thermo.D,
        W=thermo.W,
        R=thermo.R,
    )
end

@inline function omega_identity_residual(thermo, T::Real)
    return thermo.energy + thermo.pressure - Float64(T) * thermo.entropy -
        thermo.mu_p * thermo.rho_p - thermo.mu_n * thermo.rho_n
end

"""Finite-difference checks at fixed fields, kept separate from production EOS."""
function thermodynamic_consistency_report(st::GasLiquidState, T::Real, p::GasLiquidCoreParams;
    p_num::Int=96, p_max_inv_fm::Real=10.0, step::Float64=1e-5)
    dens = density_bundle(st, T, p; p_num=p_num, p_max_inv_fm=p_max_inv_fm)
    fields = field_contributions(dens, p)
    thermo = pressure_density_entropy_energy(st, T, p; p_num=p_num, p_max_inv_fm=p_max_inv_fm)
    hp = GasLiquidState(st.S, st.D, st.mu_p + step, st.mu_n)
    hn = GasLiquidState(st.S, st.D, st.mu_p, st.mu_n + step)
    pp = -_omega_fixed_fields(hp, T, p, fields.W, fields.R; p_num=p_num, p_max_inv_fm=p_max_inv_fm)
    pm = -_omega_fixed_fields(GasLiquidState(st.S, st.D, st.mu_p - step, st.mu_n), T, p, fields.W, fields.R; p_num=p_num, p_max_inv_fm=p_max_inv_fm)
    pn = -_omega_fixed_fields(hn, T, p, fields.W, fields.R; p_num=p_num, p_max_inv_fm=p_max_inv_fm)
    pnm = -_omega_fixed_fields(GasLiquidState(st.S, st.D, st.mu_p, st.mu_n - step), T, p, fields.W, fields.R; p_num=p_num, p_max_inv_fm=p_max_inv_fm)
    dP_dmu_p = (pp - pm) / (2.0 * step)
    dP_dmu_n = (pn - pnm) / (2.0 * step)
    ht = max(step, 1e-4 * max(Float64(T), 1e-3))
    dP_dT = (-_omega_fixed_fields(st, T + ht, p, fields.W, fields.R; p_num=p_num, p_max_inv_fm=p_max_inv_fm) +
        _omega_fixed_fields(st, max(T - ht, 1e-10), p, fields.W, fields.R; p_num=p_num, p_max_inv_fm=p_max_inv_fm)) / (2.0 * ht)
    return (
        rho_p_error=dP_dmu_p - dens.rho_p,
        rho_n_error=dP_dmu_n - dens.rho_n,
        entropy_error=dP_dT - thermo.entropy,
        identity_error=omega_identity_residual(thermo, T),
        finite=isfinite(dP_dmu_p) && isfinite(dP_dmu_n) && isfinite(dP_dT),
    )
end

end # module
