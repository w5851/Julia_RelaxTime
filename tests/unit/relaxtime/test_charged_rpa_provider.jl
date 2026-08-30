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
using Main.RelaxTime.ChargedRPAProvider: charged_polarization, charged_pole_residual,
    charged_mott_diagnostic
using Main.RelaxTime.PolarizationAniso: polarization_aniso, polarization_with_width

@testset "ChargedRPAProvider Phase C contract" verbose=true begin
    masses = (u=1.50, d=1.55, s=2.50)
    chemical_potentials = (u=0.20, d=-0.10, s=0.05)
    A_values = (u=-3.0, d=-3.1, s=-4.5)
    thermo = (T=0.75, Φ=0.4, Φbar=0.45, ξ=0.1)
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
        @test plus.num_s_quark == 1
        @test minus.num_s_quark == 1
        @test plus.provider == :PolarizationAniso
        @test plus.mode == :real_axis
        @test plus.retarded_convention == :external_retarded
        @test isfinite(real(plus.value)) && isfinite(imag(plus.value))
        @test isfinite(real(minus.value)) && isfinite(imag(minus.value))
        @test plus.value != minus.value
    end

    @testset "provider reuses current A/B0 algebra" begin
        spec = charged_rpa_spec(:K_plus; channel=:P)
        sample = charged_polarization(spec, k0, q, masses, chemical_potentials, thermo, A_values)
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
            thermo.ξ,
            A_values.u,
            A_values.s,
            1,
        )
        @test sample.value ≈ ComplexF64(expected_re, expected_im)
        @test sample.polarization_units == :fm_minus2

        finite_width = charged_polarization(
            spec,
            k0,
            q,
            masses,
            chemical_potentials,
            thermo,
            A_values;
            gamma_inv_fm=0.04,
            mode=:finite_width,
        )
        expected_width_re, expected_width_im = polarization_with_width(
            :P,
            k0,
            0.04,
            q,
            masses.u,
            masses.s,
            chemical_potentials.u,
            chemical_potentials.s,
            thermo.T,
            thermo.Φ,
            thermo.Φbar,
            thermo.ξ,
            A_values.u,
            A_values.s,
            1,
        )
        @test finite_width.mode == :finite_width
        @test finite_width.gamma_inv_fm == 0.04
        @test finite_width.value ≈ ComplexF64(expected_width_re, expected_width_im)
    end

    @testset "input validation and continuation boundary" begin
        spec = charged_rpa_spec(:pi_plus)
        @test_throws ArgumentError charged_polarization(spec, k0, -q, masses, chemical_potentials, thermo, A_values)
        @test_throws ArgumentError charged_polarization(spec, k0, q, masses, chemical_potentials, thermo, A_values; gamma_inv_fm=0.1)
        @test_throws ArgumentError charged_polarization(spec, k0, q, masses, chemical_potentials, thermo, A_values; mode=:retarded)
        @test_throws ArgumentError charged_polarization(spec, k0, q, masses, chemical_potentials, (T=0.0, Φ=0.4, Φbar=0.4), A_values)
        @test_throws ArgumentError charged_polarization(spec, k0, q, (u=1.5, d=1.5), chemical_potentials, thermo, A_values)
        @test_throws ArgumentError charged_polarization(spec, k0, q, masses, chemical_potentials, thermo, A_values; num_s_quark=2)
        @test_throws ArgumentError charged_polarization(spec, k0, q, masses, chemical_potentials, thermo, A_values; gamma_inv_fm=-0.1, mode=:finite_width)
    end

    @testset "pole residual and Mott branch records" begin
        spec = charged_rpa_spec(:K_plus)
        K_a = 0.35
        Pi_a = 0.20 + 0.03im
        residual = charged_pole_residual(spec, K_a, Pi_a)
        expected = 1 - 4 * K_a * Pi_a
        @test residual.residual_complex ≈ expected
        @test residual.residual_real ≈ real(expected)
        @test residual.residual_imag ≈ imag(expected)
        @test residual.residual_norm ≈ abs(expected)
        @test residual.root_search == false
        @test residual.polarization_units == :fm_minus2

        bound = charged_mott_diagnostic(spec, 3.9, masses)
        at_threshold = charged_mott_diagnostic(spec, 4.0, masses; pole_gamma_inv_fm=0.05)
        continuum = charged_mott_diagnostic(spec, 4.1, masses)
        @test bound.threshold_inv_fm == masses.u + masses.s
        @test bound.status == :bound
        @test bound.is_mott == false
        @test at_threshold.status == :at_threshold
        @test at_threshold.is_mott == true
        @test at_threshold.pole_gamma_inv_fm == 0.05
        @test continuum.status == :continuum
        @test continuum.threshold_gap_inv_fm > 0.0
        @test_throws ArgumentError charged_mott_diagnostic(spec, -0.1, masses)
        @test_throws ArgumentError charged_mott_diagnostic(spec, 4.0, masses; atol=-1e-6)
    end
end
