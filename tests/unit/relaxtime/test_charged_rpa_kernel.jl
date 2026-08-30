"""Unit tests for the Phase-B charged scalar-RPA kernel contract."""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

const _ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _RELAXTIME_PATH = joinpath(_ROOT, "src", "relaxtime", "RelaxTime.jl")
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAXTIME_PATH)
end

using Main.RelaxTime.ChargedRPAKernel
using Main.RelaxTime.EffectiveCouplings: calculate_effective_couplings_from_phi
using Main.RelaxTime.MesonInteractionKernel: build_full_kmt_interaction
using Main.RelaxTime.MesonPropagator: meson_propagator_simple

@testset "ChargedRPAKernel Phase B contract" verbose=true begin
    @testset "charged flavor and KMT mappings" begin
        pi_plus = charged_rpa_spec(:pi_plus)
        pi_minus = charged_rpa_spec(:pi_minus; channel=:pseudoscalar)
        K_plus = charged_rpa_spec(:K_plus)
        K_minus = charged_rpa_spec(:K_minus; channel=:plus)

        @test pi_plus.pair == (:u, :d)
        @test pi_minus.pair == (:d, :u)
        @test pi_plus.kernel_pair == :K12
        @test pi_minus.kernel_pair == :K12
        @test K_plus.pair == (:u, :s)
        @test K_minus.pair == (:s, :u)
        @test K_plus.kernel_pair == :K45
        @test K_minus.kernel_pair == :K45
        @test pi_plus.channel == :P
        @test K_minus.channel == :P
        @test pi_plus.normalization_source == :legacy_scalar_diagnostic
        @test pi_plus.retarded_convention == :external_retarded
    end

    @testset "P/S and normalization metadata" begin
        scalar = charged_rpa_spec(:K_plus; channel=:S)
        candidate = charged_rpa_spec(:K_plus; normalization=:matrix_scalar_candidate)

        @test scalar.channel == :S
        @test scalar.numerator_factor == 2.0
        @test scalar.denominator_factor == 4.0
        @test candidate.numerator_factor == 2.0
        @test candidate.denominator_factor == 2.0
        @test candidate.normalization_source == :neutral_matrix_candidate
        @test candidate.retarded_convention == :external_retarded
    end

    @testset "full kernel coupling selection" begin
        kernel = build_full_kmt_interaction((-0.30, -0.20, -0.10); G=1.25, K=0.40)
        K_plus = charged_rpa_spec(:K_plus)
        K_minus = charged_rpa_spec(:K_minus)
        pi = charged_rpa_spec(:pi_plus)

        @test charged_rpa_coupling(kernel, K_plus) ≈ kernel.charged_P.K45
        @test charged_rpa_coupling(K_minus, kernel) ≈ kernel.charged_P.K45
        @test charged_rpa_coupling(kernel, pi) ≈ kernel.charged_P.K12
    end

    @testset "scalar algebra is explicit and unregularized" begin
        legacy = charged_rpa_spec(:K_plus)
        candidate = charged_rpa_spec(:K_plus; normalization=:matrix_scalar_candidate)
        K_a = 0.35
        Pi_a = 0.20 + 0.03im

        @test charged_rpa_inverse(legacy, K_a, Pi_a) ≈ 1 - 4 * K_a * Pi_a
        @test charged_rpa_propagator(legacy, K_a, Pi_a) ≈ 2 * K_a / (1 - 4 * K_a * Pi_a)
        @test charged_rpa_inverse(candidate, K_a, Pi_a) ≈ 1 - 2 * K_a * Pi_a
        @test charged_rpa_propagator(candidate, K_a, Pi_a) ≈ 2 * K_a / (1 - 2 * K_a * Pi_a)
        @test charged_rpa_propagator(legacy, K_a, Pi_a) != charged_rpa_propagator(candidate, K_a, Pi_a)

        # Exact poles are exposed to the caller; no epsilon is silently added.
        @test charged_rpa_inverse(legacy, 0.25, 1.0) == 0.0
        @test_throws DomainError charged_rpa_propagator(legacy, 0.25, 1.0)
    end

    @testset "isospin-symmetric parity with legacy scalar propagator" begin
        G = 1.25
        K = 0.40
        phi_l, phi_s = -0.30, -0.10
        kernel = build_full_kmt_interaction((phi_l, phi_l, phi_s); G=G, K=K)
        legacy = calculate_effective_couplings_from_phi(G, K, phi_l, phi_s)

        cases = (
            (:pi_plus, :pi, :P, legacy.K123_plus, 0.11 + 0.02im),
            (:pi_minus, :pi, :P, legacy.K123_plus, 0.11 - 0.02im),
            (:K_plus, :K, :P, legacy.K4567_plus, 0.09 + 0.015im),
            (:K_minus, :K, :P, legacy.K4567_plus, 0.09 - 0.015im),
            (:pi_plus, :sigma_pi, :S, legacy.K123_minus, 0.07 + 0.01im),
            (:K_plus, :sigma_K, :S, legacy.K4567_minus, 0.06 + 0.01im),
        )

        for (charged_meson, legacy_meson, channel, coupling, Pi_a) in cases
            spec = charged_rpa_spec(charged_meson; channel=channel)
            @test charged_rpa_coupling(kernel, spec) ≈ coupling
            @test charged_rpa_inverse(spec, coupling, Pi_a) ≈ 1 - 4 * coupling * Pi_a
            @test charged_rpa_propagator(spec, coupling, Pi_a) ≈
                meson_propagator_simple(legacy_meson, legacy, Pi_a)
        end

        # The same denominator gives the same pole location in the parity test;
        # the strict helper exposes the exact pole instead of adding epsilon.
        spec = charged_rpa_spec(:K_plus)
        coupling = charged_rpa_coupling(kernel, spec)
        Pi_near_pole = ComplexF64((1 / (4 * coupling)) * (1 - 1e-6))
        @test charged_rpa_inverse(spec, coupling, Pi_near_pole) ≈
            1 - 4 * coupling * Pi_near_pole
        @test charged_rpa_propagator(spec, coupling, Pi_near_pole) ≈
            meson_propagator_simple(:K, legacy, Pi_near_pole)
    end

    @testset "K-plus/K-minus polarization inputs remain ordered" begin
        K_plus = charged_rpa_spec(:K_plus)
        K_minus = charged_rpa_spec(:K_minus)
        @test K_plus.pair == (:u, :s)
        @test K_minus.pair == (:s, :u)
        @test K_plus.pair != K_minus.pair
        @test K_plus.kernel_pair == K_minus.kernel_pair == :K45

        K_a = 0.35
        Pi_us = 0.20 + 0.03im
        Pi_su = 0.17 - 0.02im
        D_plus = charged_rpa_propagator(K_plus, K_a, Pi_us)
        D_minus = charged_rpa_propagator(K_minus, K_a, Pi_su)
        @test D_plus ≈ 2 * K_a / (1 - 4 * K_a * Pi_us)
        @test D_minus ≈ 2 * K_a / (1 - 4 * K_a * Pi_su)
        @test D_plus != D_minus
    end

    @testset "input validation" begin
        @test_throws ArgumentError charged_rpa_spec(:K0)
        @test_throws ArgumentError charged_rpa_spec(:pi_plus; channel=:unknown)
        @test_throws ArgumentError charged_rpa_spec(:K_plus; normalization=:unknown)
        @test_throws ArgumentError charged_rpa_spec(:K_plus; retarded_convention=Symbol(""))

        spec = charged_rpa_spec(:pi_plus)
        @test_throws ArgumentError charged_rpa_inverse(spec, Inf, 0.1 + 0im)
        @test_throws ArgumentError charged_rpa_inverse(spec, 0.2, NaN + 0im)
        @test_throws ArgumentError charged_rpa_propagator(spec, 0.2, Inf + 0im)
    end
end
