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
using Main.RelaxTime.MesonInteractionKernel: build_full_kmt_interaction

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
