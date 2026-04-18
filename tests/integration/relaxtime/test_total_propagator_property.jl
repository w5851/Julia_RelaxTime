"""
Property-based tests for TotalPropagator module struct equivalence.

Tests Property 1: Struct-NamedTuple Equivalence (TotalPropagator)
Validates: Requirements 6.2, 11.6

This test verifies that total propagator functions produce identical results
whether called with QuarkParams/ThermoParams structs or NamedTuple parameters.
"""

using Test
const INTEGRATION_VERBOSE = get(ENV, "INTEGRATION_VERBOSE", "0") in ("1", "true", "TRUE", "yes", "YES")
@inline integration_println(args...) = INTEGRATION_VERBOSE ? Base.println(args...) : nothing

using Supposition

push!(LOAD_PATH, joinpath(@__DIR__, "../../../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src/relaxtime"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src/utils"))

if !isdefined(Main, :RelaxTime)
    include("../../../src/relaxtime/RelaxTime.jl")
end

using .TotalPropagator
using .Constants_PNJL
using .EffectiveCouplings
using .OneLoopIntegrals: A
using .GaussLegendre: gauleg

include("test_utils.jl")
using .Main: QuarkParams, ThermoParams, as_namedtuple

const _ci = get(ENV, "CI", "") in ("1", "true", "TRUE", "yes", "YES")

@inline finite_float(minimum::Float64, maximum::Float64) =
    Data.Floats{Float64}(minimum=minimum, maximum=maximum, nans=false, infs=false)

@testset "TotalPropagator Property Tests" begin
    @testset "Property: total_propagator_simple Equivalence" begin
        integration_println("\n" * "="^70)
        integration_println("Property Test: total_propagator_simple Struct-NamedTuple Equivalence")
        integration_println("Testing with random parameters")
        integration_println("="^70)

        RTOL = 1e-12
        ATOL = 1e-14
        nodes_p, weights_p = gauleg(0.0, 20.0, 64)

        @check max_examples=(_ci ? 10 : 20) function property_total_propagator_simple_equivalence(
            m_u = finite_float(0.5, 2.0),
            m_s = finite_float(2.0, 5.0),
            μ_u = finite_float(0.0, 0.5),
            μ_s = finite_float(0.0, 0.5),
            T = finite_float(0.1, 0.3),
            Φ = finite_float(0.1, 0.9),
            Φbar = finite_float(0.1, 0.9),
            k0 = finite_float(0.1, 2.0),
            k_norm = finite_float(0.0, 1.5),
        )
            A_u = A(m_u, μ_u, T, Φ, Φbar, nodes_p, weights_p)
            A_s = A(m_s, μ_s, T, Φ, Φbar, nodes_p, weights_p)
            G_u = EffectiveCouplings.calculate_G_from_A(A_u, m_u)
            G_s = EffectiveCouplings.calculate_G_from_A(A_s, m_s)
            K_coeffs = EffectiveCouplings.calculate_effective_couplings(Constants_PNJL.G_fm2, Constants_PNJL.K_fm5, G_u, G_s)

            q_nt = (m=(u=m_u, d=m_u, s=m_s), μ=(u=μ_u, d=μ_u, s=μ_s), A=(u=A_u, d=A_u, s=A_s))
            t_nt = (T=T, Φ=Φ, Φbar=Φbar, ξ=0.0)

            q_struct = QuarkParams((u=m_u, d=m_u, s=m_s), (u=μ_u, d=μ_u, s=μ_s))
            t_struct = ThermoParams(T, Φ, Φbar, 0.0)
            q_struct_with_A = merge(as_namedtuple(q_struct), (A=(u=A_u, d=A_u, s=A_s),))

            D_struct = total_propagator_simple(:uu_to_uu, :t, [:pi, :K], k0, k_norm, q_struct_with_A, t_struct, K_coeffs)
            D_nt = total_propagator_simple(:uu_to_uu, :t, [:pi, :K], k0, k_norm, q_nt, t_nt, K_coeffs)

            isapprox(real(D_struct), real(D_nt), rtol=RTOL, atol=ATOL) || error("Real part mismatch")
            isapprox(imag(D_struct), imag(D_nt), rtol=RTOL, atol=ATOL) || error("Imag part mismatch")
            isfinite(D_struct) || error("Propagator must be finite (struct): $D_struct")
            isfinite(D_nt) || error("Propagator must be finite (nt): $D_nt")
        end

        integration_println("✓ Property test passed: total_propagator_simple struct-NamedTuple equivalence verified")
    end

    @testset "Property: total_propagator_mixed Equivalence" begin
        integration_println("\n" * "="^70)
        integration_println("Property Test: total_propagator_mixed Struct-NamedTuple Equivalence")
        integration_println("Testing with random parameters")
        integration_println("="^70)

        RTOL = 1e-12
        ATOL = 1e-14
        nodes_p, weights_p = gauleg(0.0, 20.0, 64)

        @check max_examples=(_ci ? 10 : 20) function property_total_propagator_mixed_equivalence(
            m_u = finite_float(0.5, 2.0),
            m_s = finite_float(2.0, 5.0),
            μ_u = finite_float(0.0, 0.5),
            μ_s = finite_float(0.0, 0.5),
            T = finite_float(0.1, 0.3),
            Φ = finite_float(0.1, 0.9),
            Φbar = finite_float(0.1, 0.9),
            k0 = finite_float(0.1, 2.0),
            k_norm = finite_float(0.0, 1.5),
        )
            A_u = A(m_u, μ_u, T, Φ, Φbar, nodes_p, weights_p)
            A_s = A(m_s, μ_s, T, Φ, Φbar, nodes_p, weights_p)
            G_u = EffectiveCouplings.calculate_G_from_A(A_u, m_u)
            G_s = EffectiveCouplings.calculate_G_from_A(A_s, m_s)
            K_coeffs = EffectiveCouplings.calculate_effective_couplings(Constants_PNJL.G_fm2, Constants_PNJL.K_fm5, G_u, G_s)

            q_struct = QuarkParams((u=m_u, d=m_u, s=m_s), (u=μ_u, d=μ_u, s=μ_s))
            t_struct = ThermoParams(T, Φ, Φbar, 0.0)
            q_struct_with_A = merge(as_namedtuple(q_struct), (A=(u=A_u, d=A_u, s=A_s),))

            q_nt = (m=(u=m_u, d=m_u, s=m_s), μ=(u=μ_u, d=μ_u, s=μ_s), A=(u=A_u, d=A_u, s=A_s))
            t_nt = (T=T, Φ=Φ, Φbar=Φbar, ξ=0.0)

            D_struct = total_propagator_mixed(:uu_to_uu, :t, :P, k0, k_norm, q_struct_with_A, t_struct, K_coeffs)
            D_nt = total_propagator_mixed(:uu_to_uu, :t, :P, k0, k_norm, q_nt, t_nt, K_coeffs)
            isapprox(real(D_struct), real(D_nt), rtol=RTOL, atol=ATOL) || error("Real part mismatch (P channel)")
            isapprox(imag(D_struct), imag(D_nt), rtol=RTOL, atol=ATOL) || error("Imag part mismatch (P channel)")
            isfinite(D_struct) || error("Propagator must be finite (struct)")
            isfinite(D_nt) || error("Propagator must be finite (nt)")

            D_struct_S = total_propagator_mixed(:uu_to_uu, :t, :S, k0, k_norm, q_struct_with_A, t_struct, K_coeffs)
            D_nt_S = total_propagator_mixed(:uu_to_uu, :t, :S, k0, k_norm, q_nt, t_nt, K_coeffs)
            isapprox(real(D_struct_S), real(D_nt_S), rtol=RTOL, atol=ATOL) || error("Real part mismatch (S channel)")
            isapprox(imag(D_struct_S), imag(D_nt_S), rtol=RTOL, atol=ATOL) || error("Imag part mismatch (S channel)")
        end

        integration_println("✓ Property test passed: total_propagator_mixed struct-NamedTuple equivalence verified")
    end

    @testset "Property: calculate_cms_momentum Equivalence" begin
        integration_println("\n" * "="^70)
        integration_println("Property Test: calculate_cms_momentum Struct-NamedTuple Equivalence")
        integration_println("Testing with random kinematic variables")
        integration_println("="^70)

        RTOL = 1e-12
        ATOL = 1e-14

        @check max_examples=(_ci ? 15 : 30) function property_calculate_cms_momentum_equivalence(
            m_u = finite_float(0.5, 2.0),
            m_s = finite_float(2.0, 5.0),
            μ_u = finite_float(0.0, 0.5),
            μ_s = finite_float(0.0, 0.5),
            s = finite_float(5.0, 15.0),
            t = finite_float(-1.5, -0.2),
        )
            q_struct = QuarkParams((u=m_u, d=m_u, s=m_s), (u=μ_u, d=μ_u, s=μ_s))
            q_nt = (m=(u=m_u, d=m_u, s=m_s), μ=(u=μ_u, d=μ_u, s=μ_s))

            result_struct = calculate_cms_momentum(:uu_to_uu, s, t, :t, q_struct)
            result_nt = calculate_cms_momentum(:uu_to_uu, s, t, :t, q_nt)
            isapprox(result_struct.k0, result_nt.k0, rtol=RTOL, atol=ATOL) || error("k0 mismatch")
            isapprox(result_struct.k, result_nt.k, rtol=RTOL, atol=ATOL) || error("k mismatch")
            result_struct.k0 >= 0.0 || error("k0 must be non-negative (struct)")
            result_struct.k >= 0.0 || error("k must be non-negative (struct)")
            result_nt.k0 >= 0.0 || error("k0 must be non-negative (nt)")
            result_nt.k >= 0.0 || error("k must be non-negative (nt)")
            isfinite(result_struct.k0) || error("k0 must be finite (struct)")
            isfinite(result_struct.k) || error("k must be finite (struct)")
            isfinite(result_nt.k0) || error("k0 must be finite (nt)")
            isfinite(result_nt.k) || error("k must be finite (nt)")
        end

        integration_println("✓ Property test passed: calculate_cms_momentum struct-NamedTuple equivalence verified")
    end

    @testset "Property: Multiple Channels and Processes" begin
        integration_println("\n" * "="^70)
        integration_println("Property Test: Multiple Channels and Processes")
        integration_println("Testing struct-NamedTuple equivalence across representative valid configurations")
        integration_println("="^70)

        RTOL = 1e-12
        ATOL = 1e-14
        nodes_p, weights_p = gauleg(0.0, 20.0, 64)

        simple_configs = [
            (:uu_to_uu, :t, [:pi, :K]),
            (:uu_to_uu, :u, [:pi, :K]),
            (:ud_to_ud, :t, [:pi]),
            (:ud_to_ud, :s, [:pi]),
            (:us_to_us, :t, [:K]),
            (:uubar_to_uubar, :s, [:pi, :K]),
            (:uubar_to_ddbar, :s, [:pi]),
        ]
        mixed_configs = [
            (:uu_to_uu, :t, :P),
            (:uu_to_uu, :t, :S),
            (:uubar_to_uubar, :s, :P),
            (:uubar_to_uubar, :s, :S),
            (:ubarsbar_to_ubarsbar, :s, :P),
            (:ubarsbar_to_ubarsbar, :s, :S),
        ]

        representative_points = [
            (m_u=1.2, m_s=3.3, μ_u=0.25, μ_s=0.28, T=0.17, Φ=0.42, Φbar=0.39, k0=0.8, k_norm=0.45),
            (m_u=1.45, m_s=3.6, μ_u=0.31, μ_s=0.34, T=0.20, Φ=0.55, Φbar=0.52, k0=1.1, k_norm=0.70),
            (m_u=1.7, m_s=3.9, μ_u=0.38, μ_s=0.36, T=0.23, Φ=0.63, Φbar=0.60, k0=1.35, k_norm=0.92),
        ]

        for point in representative_points
            A_u = A(point.m_u, point.μ_u, point.T, point.Φ, point.Φbar, nodes_p, weights_p)
            A_s = A(point.m_s, point.μ_s, point.T, point.Φ, point.Φbar, nodes_p, weights_p)
            G_u = EffectiveCouplings.calculate_G_from_A(A_u, point.m_u)
            G_s = EffectiveCouplings.calculate_G_from_A(A_s, point.m_s)
            K_coeffs = EffectiveCouplings.calculate_effective_couplings(Constants_PNJL.G_fm2, Constants_PNJL.K_fm5, G_u, G_s)

            q_struct = QuarkParams((u=point.m_u, d=point.m_u, s=point.m_s), (u=point.μ_u, d=point.μ_u, s=point.μ_s))
            t_struct = ThermoParams(point.T, point.Φ, point.Φbar, 0.0)
            q_struct_with_A = merge(as_namedtuple(q_struct), (A=(u=A_u, d=A_u, s=A_s),))
            q_nt = (m=(u=point.m_u, d=point.m_u, s=point.m_s), μ=(u=point.μ_u, d=point.μ_u, s=point.μ_s), A=(u=A_u, d=A_u, s=A_s))
            t_nt = (T=point.T, Φ=point.Φ, Φbar=point.Φbar, ξ=0.0)

            for (process, channel, meson_list) in simple_configs
                D_struct = total_propagator_simple(process, channel, meson_list, point.k0, point.k_norm, q_struct_with_A, t_struct, K_coeffs)
                D_nt = total_propagator_simple(process, channel, meson_list, point.k0, point.k_norm, q_nt, t_nt, K_coeffs)
                @test isfinite(D_struct)
                @test isfinite(D_nt)
                @test isapprox(real(D_struct), real(D_nt), rtol=RTOL, atol=ATOL)
                @test isapprox(imag(D_struct), imag(D_nt), rtol=RTOL, atol=ATOL)
            end

            for (process, channel, meson_channel) in mixed_configs
                D_struct = total_propagator_mixed(process, channel, meson_channel, point.k0, point.k_norm, q_struct_with_A, t_struct, K_coeffs)
                D_nt = total_propagator_mixed(process, channel, meson_channel, point.k0, point.k_norm, q_nt, t_nt, K_coeffs)
                @test isfinite(D_struct)
                @test isfinite(D_nt)
                @test isapprox(real(D_struct), real(D_nt), rtol=RTOL, atol=ATOL)
                @test isapprox(imag(D_struct), imag(D_nt), rtol=RTOL, atol=ATOL)
            end
        end

        integration_println("✓ Property test passed: representative multi-config propagator equivalence verified")
    end

    @testset "Property: Normalization Helpers" begin
        integration_println("\n" * "="^70)
        integration_println("Property Test: Normalization Helper Functions")
        integration_println("Testing ParameterAdapters normalization helpers")
        integration_println("="^70)

        representative_points = [
            (m_u=0.8, m_s=2.8, μ_u=0.0, μ_s=0.1, T=0.08, Φ=0.1, Φbar=0.2),
            (m_u=1.3, m_s=3.7, μ_u=0.3, μ_s=0.4, T=0.18, Φ=0.5, Φbar=0.55),
            (m_u=1.9, m_s=4.8, μ_u=0.8, μ_s=0.9, T=0.28, Φ=0.9, Φbar=0.85),
        ]

        isdefined(Main, :ParameterAdapters) || error("ParameterAdapters not defined")

        for point in representative_points
            q_struct = QuarkParams((u=point.m_u, d=point.m_u, s=point.m_s), (u=point.μ_u, d=point.μ_u, s=point.μ_s))
            t_struct = ThermoParams(point.T, point.Φ, point.Φbar, 0.0)

            q_nt = Main.ParameterAdapters.normalize_quark_input(q_struct)
            t_nt = Main.ParameterAdapters.normalize_thermo_input(t_struct)

            @test q_nt.m.u == point.m_u
            @test q_nt.m.s == point.m_s
            @test q_nt.μ.u == point.μ_u
            @test q_nt.μ.s == point.μ_s
            @test t_nt.T == point.T
            @test t_nt.Φ == point.Φ
            @test t_nt.Φbar == point.Φbar

            q_struct_2 = QuarkParams(q_nt)
            q_nt_2 = as_namedtuple(q_struct_2)
            @test q_nt == q_nt_2

            t_struct_2 = ThermoParams(t_nt)
            t_nt_2 = as_namedtuple(t_struct_2)
            @test t_nt == t_nt_2
        end

        integration_println("✓ Property test passed: Normalization helpers work correctly")
    end
end

integration_println("\n" * "="^70)
integration_println("TotalPropagator Property Tests Complete!")
integration_println("="^70)
integration_println("\nAll enabled property tests passed:")
integration_println("  1. total_propagator_simple struct-NamedTuple equivalence")
integration_println("  2. total_propagator_mixed struct-NamedTuple equivalence (P and S channels)")
integration_println("  3. calculate_cms_momentum struct-NamedTuple equivalence")
integration_println("  4. Representative multiple processes and channels equivalence")
integration_println("  5. Normalization helpers verified")
integration_println("  6. Physical constraints maintained (non-negative momenta, finite values)")
integration_println("="^70)
