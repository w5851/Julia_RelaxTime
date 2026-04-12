"""
Integration tests for TotalPropagator edge cases.

Focus:
- channel/process compatibility
- struct/NamedTuple parity on valid calls
- explicit contract checks for invalid channel and missing A field
"""

using Test

push!(LOAD_PATH, joinpath(@__DIR__, "../../../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src/relaxtime"))

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

@testset "TotalPropagator Edge Cases" begin
    m_u = 1.52
    m_s = 3.04
    μ_u = 0.3
    μ_s = 0.3
    T = 0.15
    Φ = 0.5
    Φbar = 0.5

    nodes_p, weights_p = gauleg(0.0, 20.0, 64)
    A_u = A(m_u, μ_u, T, Φ, Φbar, nodes_p, weights_p)
    A_s = A(m_s, μ_s, T, Φ, Φbar, nodes_p, weights_p)
    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)
    K_coeffs = calculate_effective_couplings(Constants_PNJL.G_fm2, Constants_PNJL.K_fm5, G_u, G_s)

    q_struct = QuarkParams((u=m_u, d=m_u, s=m_s), (u=μ_u, d=μ_u, s=μ_s))
    t_struct = ThermoParams(T, Φ, Φbar, 0.0)
    q_struct_with_A = merge(as_namedtuple(q_struct), (A=(u=A_u, d=A_u, s=A_s),))

    q_nt_with_A = (m=(u=m_u, d=m_u, s=m_s), μ=(u=μ_u, d=μ_u, s=μ_s), A=(u=A_u, d=A_u, s=A_s))
    q_nt_without_A = (m=(u=m_u, d=m_u, s=m_s), μ=(u=μ_u, d=μ_u, s=μ_s))
    t_nt = (T=T, Φ=Φ, Φbar=Φbar, ξ=0.0)

    @testset "calculate_all_propagators valid channel sets" begin
        qq = calculate_all_propagators(:uu_to_uu, 1.0, 1.0, q_struct_with_A, t_struct, K_coeffs)
        qq_nt = calculate_all_propagators(:uu_to_uu, 1.0, 1.0, q_nt_with_A, t_nt, K_coeffs)
        @test haskey(qq, :t)
        @test haskey(qq, :u)
        @test !haskey(qq, :s)
        @test isapprox(qq.t, qq_nt.t, rtol=1e-12, atol=1e-14)
        @test isapprox(qq.u, qq_nt.u, rtol=1e-12, atol=1e-14)

        qqbar = calculate_all_propagators(:uubar_to_uubar, 1.0, 1.0, q_struct_with_A, t_struct, K_coeffs)
        qqbar_nt = calculate_all_propagators(:uubar_to_uubar, 1.0, 1.0, q_nt_with_A, t_nt, K_coeffs)
        @test haskey(qqbar, :t)
        @test haskey(qqbar, :s)
        @test !haskey(qqbar, :u)
        @test isapprox(qqbar.t, qqbar_nt.t, rtol=1e-12, atol=1e-14)
        @test isapprox(qqbar.s, qqbar_nt.s, rtol=1e-12, atol=1e-14)
    end

    @testset "total_propagator_auto channel guards and parity" begin
        valid_struct = total_propagator_auto(:uu_to_uu, :t, 1.0, 1.0, q_struct_with_A, t_struct, K_coeffs)
        valid_nt = total_propagator_auto(:uu_to_uu, :t, 1.0, 1.0, q_nt_with_A, t_nt, K_coeffs)
        @test isfinite(valid_struct)
        @test isfinite(valid_nt)
        @test isapprox(valid_struct, valid_nt, rtol=1e-12, atol=1e-14)

        @test_throws ErrorException total_propagator_auto(:uu_to_uu, :s, 1.0, 1.0, q_struct_with_A, t_struct, K_coeffs)
        @test_throws ErrorException total_propagator_auto(:uubar_to_uubar, :u, 1.0, 1.0, q_struct_with_A, t_struct, K_coeffs)
    end

    @testset "calculate_cms_momentum return shape" begin
        cms_struct = calculate_cms_momentum(:uu_to_uu, 31.0, -5.0, :t, q_struct)
        cms_nt = calculate_cms_momentum(:uu_to_uu, 31.0, -5.0, :t, q_nt_without_A)
        @test haskey(cms_struct, :k0)
        @test haskey(cms_struct, :k)
        @test cms_struct.k0 >= 0.0
        @test cms_struct.k >= 0.0
        @test isapprox(cms_struct.k0, cms_nt.k0, rtol=1e-12, atol=1e-14)
        @test isapprox(cms_struct.k, cms_nt.k, rtol=1e-12, atol=1e-14)
    end

    @testset "calculate_all_propagators_by_channel requires A" begin
        by_channel_struct = calculate_all_propagators_by_channel(:uu_to_uu, 1.0, 1.0, q_struct_with_A, t_struct, K_coeffs)
        by_channel_nt = calculate_all_propagators_by_channel(:uu_to_uu, 1.0, 1.0, q_nt_with_A, t_nt, K_coeffs)

        @test haskey(by_channel_struct, :t_S)
        @test haskey(by_channel_struct, :t_P)
        @test haskey(by_channel_struct, :u_S)
        @test haskey(by_channel_struct, :u_P)
        @test isapprox(by_channel_struct.t_S, by_channel_nt.t_S, rtol=1e-12, atol=1e-14)
        @test isapprox(by_channel_struct.t_P, by_channel_nt.t_P, rtol=1e-12, atol=1e-14)
        @test isapprox(by_channel_struct.u_S, by_channel_nt.u_S, rtol=1e-12, atol=1e-14)
        @test isapprox(by_channel_struct.u_P, by_channel_nt.u_P, rtol=1e-12, atol=1e-14)

        @test_throws FieldError calculate_all_propagators_by_channel(:uu_to_uu, 1.0, 1.0, q_nt_without_A, t_nt, K_coeffs)
    end
end
