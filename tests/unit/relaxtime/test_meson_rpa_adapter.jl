"""Unit tests for the diagnostic PolarizationAniso -> neutral RPA adapter."""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

const _ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _RELAXTIME_PATH = joinpath(_ROOT, "src", "relaxtime", "RelaxTime.jl")
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAXTIME_PATH)
end

using Main.RelaxTime.MesonInteractionKernel: build_full_kmt_interaction
using Main.RelaxTime.MesonRPA: neutral_polarization_matrix, neutral_rpa_inverse_matrix,
    neutral_rpa_propagator, neutral_rpa_determinant
using Main.RelaxTime.MesonRPAAdapter: neutral_flavor_bubbles, neutral_rpa_from_quark_params
using Main.RelaxTime.PolarizationAniso: polarization_aniso, polarization_with_width

const _Q = (
    m=(u=0.36, d=0.38, s=0.55),
    μ=(u=0.08, d=0.04, s=-0.02),
    A=(u=0.17, d=0.16, s=0.14),
)
const _THERMO = (T=0.18, Φ=0.45, Φbar=0.45, ξ=0.0)

@inline function _direct_bubble(channel, flavor, k0, k_norm, num_s)
    mass = _Q.m[flavor]
    mu = _Q.μ[flavor]
    A_value = _Q.A[flavor]
    re, im = polarization_aniso(
        channel, k0, k_norm, mass, mass, mu, mu,
        _THERMO.T, _THERMO.Φ, _THERMO.Φbar, _THERMO.ξ,
        A_value, A_value, num_s,
    )
    return ComplexF64(re, im)
end

@testset "MesonRPAAdapter diagnostic bridge" verbose=true begin
    k0 = 0.42
    k_norm = 0.17

    @testset "same-flavor bubbles preserve PolarizationAniso values" begin
        result = neutral_flavor_bubbles(:P, k0, k_norm, _Q, _THERMO; ensure_A=false)

        @test (result.channel, result.k0_inv_fm, result.k_norm_inv_fm) == (:P, k0, k_norm)
        @test result.gamma_inv_fm == 0.0
        @test result.with_width == false
        @test result.a_auto_built == false
        @test result.num_s_quark == (u=0, d=0, s=0)
        for flavor in (:u, :d, :s)
            @test getproperty(result, flavor) ≈ _direct_bubble(:P, flavor, k0, k_norm, 0)
            @test isfinite(real(getproperty(result, flavor)))
            @test isfinite(imag(getproperty(result, flavor)))
        end
        @test result.quark_params.A == _Q.A
        @test result.thermo_params == _THERMO
    end

    @testset "num_s_quark and width options are forwarded" begin
        num_s = (u=0, d=1, s=2)
        result = neutral_flavor_bubbles(
            :pseudoscalar, k0, k_norm, _Q, _THERMO;
            ensure_A=false,
            num_s_quark=num_s,
        )
        @test result.channel == :P
        @test result.num_s_quark == num_s
        @test result.d ≈ _direct_bubble(:P, :d, k0, k_norm, 1)

        gamma = 0.06
        width_result = neutral_flavor_bubbles(
            :S, k0, k_norm, _Q, _THERMO;
            ensure_A=false,
            gamma_inv_fm=gamma,
            with_width=true,
        )
        re, im = polarization_with_width(
            :S, k0, gamma, k_norm,
            _Q.m.u, _Q.m.u, _Q.μ.u, _Q.μ.u,
            _THERMO.T, _THERMO.Φ, _THERMO.Φbar, _THERMO.ξ,
            _Q.A.u, _Q.A.u, 0,
        )
        @test width_result.u ≈ ComplexF64(re, im)
        @test width_result.gamma_inv_fm == gamma
        @test width_result.with_width
    end

    @testset "adapter composes with neutral RPA algebra" begin
        kernel = build_full_kmt_interaction((-0.30, -0.20, -0.10); G=1.25, K=0.40)
        result = neutral_rpa_from_quark_params(
            kernel, k0, k_norm, _Q, _THERMO;
            ensure_A=false,
        )
        pi = neutral_polarization_matrix(result.bubbles)
        @test result.polarization == pi
        @test result.inverse_matrix == neutral_rpa_inverse_matrix(kernel, pi; channel=:P)
        @test result.propagator == neutral_rpa_propagator(kernel, pi; channel=:P)
        @test result.determinant == neutral_rpa_determinant(kernel, pi; channel=:P)
        @test result.channel == :P
        @test all(isfinite, real.(result.propagator))
        @test all(isfinite, imag.(result.propagator))
    end

    @testset "automatic A and input validation" begin
        q_without_A = (m=_Q.m, μ=_Q.μ)
        auto = neutral_flavor_bubbles(
            :P, 0.35, 0.12, q_without_A, _THERMO;
            a_p_nodes=4,
            a_p_max=4.0,
            a_cos_nodes=4,
            warn_on_auto_A=false,
        )
        @test auto.a_auto_built
        @test all(flavor -> hasproperty(auto.quark_params.A, flavor), (:u, :d, :s))

        @test_throws ArgumentError neutral_flavor_bubbles(:X, k0, k_norm, _Q, _THERMO; ensure_A=false)
        @test_throws ArgumentError neutral_flavor_bubbles(:P, k0, -0.1, _Q, _THERMO; ensure_A=false)
        @test_throws ArgumentError neutral_flavor_bubbles(:P, k0, k_norm, _Q, _THERMO; ensure_A=false, gamma_inv_fm=0.1)
        @test_throws ArgumentError neutral_flavor_bubbles(:P, k0, k_norm, _Q, _THERMO; ensure_A=false, num_s_quark=(u=0, d=1.0, s=0))
        @test_throws ArgumentError neutral_flavor_bubbles(:P, k0, k_norm, (m=_Q.m, μ=_Q.μ), _THERMO; ensure_A=false)
    end
end
