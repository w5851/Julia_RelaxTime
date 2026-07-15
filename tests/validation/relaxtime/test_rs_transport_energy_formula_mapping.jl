using Test

const _RS_TRANSPORT_RELAXTIME_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RS_TRANSPORT_RELAXTIME_PATH)
end
const RS_TRANSPORT_TC = Main.TransportCoefficients

@testset "RS transport energy formula mapping" begin
    @testset "muB=0 bulk kernel matches Mykhaylova-Sasaki Eq. (10) bracket" begin
        p = 0.8
        m = 0.35
        T = 0.16
        c_s_sq = 0.24
        dM_dT = -0.12
        E_kin = sqrt(p^2 + m^2)

        B = RS_TRANSPORT_TC._bulk_isentropic_B(
            p,
            m,
            0.0,
            T,
            c_s_sq,
            0.0,
            dM_dT,
            0.0,
            false,
            E_kin,
        )

        # For Π=m², T² ∂Π/∂T² = T*m*dM/dT = T*E*dE/dT.
        mykhaylova_sasaki_bracket = c_s_sq * (E_kin^2 - T * m * dM_dT) - p^2 / 3.0
        @test isapprox(B, -3.0 * mykhaylova_sasaki_bracket; rtol=1e-12, atol=0.0)
        @test isapprox(B^2, 9.0 * mykhaylova_sasaki_bracket^2; rtol=1e-12, atol=0.0)
    end

    @testset "kappa_BQ matches Das et al. Eq. (55) energy roles" begin
        quark_params = (m=(u=0.3, d=0.3, s=0.5), μ=(u=0.1, d=0.1, s=0.1))
        thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.3)
        tau = (u=1.0, d=1.1, s=1.2, ubar=1.3, dbar=1.4, sbar=1.5)
        densities = (u=0.12, d=0.09, s=0.03, ubar=0.02, dbar=0.01, sbar=0.005)
        pressure = 0.08
        energy_density = 0.60
        p = 1.1
        p_weight = 0.6
        cosθ = 0.4
        cos_weight = 1.7
        E_kin = 2.0
        E_dist = 5.0
        occupancy = E_dist / 10.0
        degeneracy = RS_TRANSPORT_TC.degeneracy_default()
        config = RS_TRANSPORT_TC.TransportIntegrationConfig(
            p_nodes=1,
            p_max=2.0,
            p_grid=[p],
            p_w=[p_weight],
            cos_nodes=1,
            cos_grid=[cosθ],
            cos_w=[cos_weight],
        )
        provider = (
            energy_from_p=(p::Float64, m::Float64) -> E_kin,
            energy_from_p_aniso=(p::Float64, m::Float64, ξ::Float64, c::Float64) -> E_dist,
            quark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> E / 10.0,
            antiquark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> E / 10.0,
            quark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> occupancy,
            antiquark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> occupancy,
            prefer_energy_aniso=true,
        )

        actual = RS_TRANSPORT_TC.diffusion_coefficient(
            quark_params,
            thermo_params;
            tau=tau,
            charge_left=:B,
            charge_right=:Q,
            densities=densities,
            pressure=pressure,
            energy=energy_density,
            degeneracy=degeneracy,
            provider=provider,
            config=config,
        )

        charge_densities = RS_TRANSPORT_TC.conserved_charge_densities(densities)
        enthalpy = pressure + energy_density
        species = (:u, :d, :s, :ubar, :dbar, :sbar)
        species_sum = 0.0
        for sp in species
            projection_B = RS_TRANSPORT_TC.conserved_charge_value(sp, :B) - charge_densities.B * E_kin / enthalpy
            projection_Q = RS_TRANSPORT_TC.conserved_charge_value(sp, :Q) - charge_densities.Q * E_kin / enthalpy
            species_sum += p^4 / E_kin^2 * degeneracy * getproperty(tau, sp) * occupancy * projection_B * projection_Q
        end
        measure = p_weight * cos_weight / (4.0 * π^2)
        expected = measure * species_sum / (3.0 * thermo_params.T^2)

        @test isapprox(actual, expected; rtol=1e-12, atol=0.0)
    end
end
