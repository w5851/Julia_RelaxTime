"""Unit tests for the Phase-C charged polarization provider contract."""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

const _ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _RELAXTIME_PATH = joinpath(_ROOT, "src", "relaxtime", "RelaxTime.jl")
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAXTIME_PATH)
end

using Main.RelaxTime.ChargedRPAKernel: charged_rpa_spec
using Main.RelaxTime.ChargedRPAProvider: charged_polarization
using Main.RelaxTime.PolarizationAniso: polarization_aniso

@testset "ChargedRPAProvider Phase C contract" verbose=true begin
    masses = (u=1.50, d=1.55, s=2.50)
    chemical_potentials = (u=0.20, d=-0.10, s=0.05)
    A_values = (u=-3.0, d=-3.1, s=-4.5)
    thermo = (T=0.75, Φ=0.4, Φbar=0.45, ξ=0.0)
    thermo_aniso = merge(thermo, (ξ=0.1,))
    k0, q = 0.80, 0.35

    @testset "ordered K-plus/K-minus inputs" begin
        plus = charged_polarization(
            charged_rpa_spec(:K_plus),
            k0,
            q,
            masses,
            chemical_potentials,
            thermo,
            A_values,
        )
        minus = charged_polarization(
            charged_rpa_spec(:K_minus),
            k0,
            q,
            masses,
            chemical_potentials,
            thermo,
            A_values,
        )

        @test plus.pair == (:u, :s)
        @test minus.pair == (:s, :u)
        @test plus.kernel_pair == minus.kernel_pair == :K45
        @test plus.m1_inv_fm == masses.u
        @test plus.m2_inv_fm == masses.s
        @test minus.m1_inv_fm == masses.s
        @test minus.m2_inv_fm == masses.u
        @test plus.μ1_inv_fm == chemical_potentials.u
        @test plus.μ2_inv_fm == chemical_potentials.s
        @test minus.μ1_inv_fm == chemical_potentials.s
        @test minus.μ2_inv_fm == chemical_potentials.u
        @test plus.A1_inv_fm2 == A_values.u
        @test plus.A2_inv_fm2 == A_values.s
        @test minus.A1_inv_fm2 == A_values.s
        @test minus.A2_inv_fm2 == A_values.u
        @test plus.num_s_quark == 0
        @test minus.num_s_quark == 0
        @test plus.provider == :OneLoopIntegralsRetarded
        @test plus.prescription == :ordered_retarded
        @test plus.analytic_scope == :upper_half_plane_probe
        @test plus.retarded_convention == :retarded_e_minus_iwt
        @test plus.eta_inv_fm == 1.0e-3
        @test plus.energy_nodes == 128
        @test isfinite(real(plus.value)) && isfinite(imag(plus.value))
        @test isfinite(real(minus.value)) && isfinite(imag(minus.value))
        @test plus.value != minus.value
    end

    @testset "strict ordered relation and legacy oracles" begin
        spec = charged_rpa_spec(:K_plus; channel=:P)
        sample = charged_polarization(
            spec, k0, q, masses, chemical_potentials, thermo, A_values;
            eta_inv_fm=0.003,
            energy_nodes=128,
        )
        conjugate_reverse = charged_polarization(
            charged_rpa_spec(:K_minus; channel=:P),
            -k0,
            q,
            masses,
            chemical_potentials,
            thermo,
            A_values;
            eta_inv_fm=0.003,
            energy_nodes=128,
        )
        @test sample.value ≈ conj(conjugate_reverse.value) rtol=1e-12 atol=1e-12
        @test sample.polarization_units == :fm_minus2

        ordered_legacy = charged_polarization(
            spec,
            k0,
            q,
            masses,
            chemical_potentials,
            thermo_aniso,
            A_values;
            prescription=:ordered_legacy_B0,
        )
        expected_re, expected_im = polarization_aniso(
            :P,
            k0,
            q,
            masses.u,
            masses.s,
            chemical_potentials.u,
            chemical_potentials.s,
            thermo.T,
            thermo.Φ,
            thermo.Φbar,
            thermo_aniso.ξ,
            A_values.u,
            A_values.s,
            0,
        )
        @test ordered_legacy.value ≈ ComplexF64(expected_re, expected_im)
        @test ordered_legacy.provider == :PolarizationAniso
        @test ordered_legacy.prescription == :ordered_legacy_B0
        @test ordered_legacy.analytic_scope == :real_axis_legacy
        @test ordered_legacy.eta_inv_fm == 0.0
        @test ordered_legacy.energy_nodes == 0

        legacy = charged_polarization(
            spec,
            k0,
            q,
            masses,
            chemical_potentials,
            thermo_aniso,
            A_values;
            prescription=:legacy_symmetrized_B0,
        )
        expected_legacy_re, expected_legacy_im = polarization_aniso(
            :P,
            k0,
            q,
            masses.u,
            masses.s,
            chemical_potentials.u,
            chemical_potentials.s,
            thermo.T,
            thermo.Φ,
            thermo.Φbar,
            thermo_aniso.ξ,
            A_values.u,
            A_values.s,
            1,
        )
        @test legacy.prescription == :legacy_symmetrized_B0
        @test legacy.num_s_quark == 1
        @test legacy.value ≈ ComplexF64(expected_legacy_re, expected_legacy_im)
        @test legacy.value != ordered_legacy.value

        pion_legacy = charged_polarization(
            charged_rpa_spec(:pi_plus),
            k0,
            q,
            masses,
            chemical_potentials,
            thermo_aniso,
            A_values;
            prescription=:legacy_symmetrized_B0,
        )
        @test pion_legacy.num_s_quark == 0
    end

    @testset "input validation and prescription boundary" begin
        spec = charged_rpa_spec(:pi_plus)
        @test_throws ArgumentError charged_polarization(spec, k0, -q, masses, chemical_potentials, thermo, A_values)
        @test_throws ArgumentError charged_polarization(spec, k0, q, masses, chemical_potentials, thermo, A_values; prescription=:unknown)
        @test_throws ArgumentError charged_polarization(spec, k0, q, masses, chemical_potentials, (T=0.0, Φ=0.4, Φbar=0.4), A_values)
        @test_throws ArgumentError charged_polarization(spec, k0, q, (u=1.5, d=1.5), chemical_potentials, thermo, A_values)
        @test_throws ArgumentError charged_polarization(spec, k0, q, masses, chemical_potentials, thermo_aniso, A_values)
        @test_throws ArgumentError charged_polarization(spec, k0, q, masses, chemical_potentials, thermo, A_values; eta_inv_fm=0.0)
        @test_throws ArgumentError charged_polarization(spec, k0, q, masses, chemical_potentials, thermo, A_values; energy_nodes=3)
    end
end
