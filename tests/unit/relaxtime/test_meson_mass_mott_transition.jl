"""
MesonMass 与 MottTransition 基础测试（轻量验证）。
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using Test

include("../../../src/Constants_PNJL.jl")
include("../../../src/relaxtime/EffectiveCouplings.jl")
include("../../../src/relaxtime/MesonMass.jl")
include("../../../src/relaxtime/MottTransition.jl")

using .Constants_PNJL: G_fm2, K_fm5
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .MesonMass: ensure_quark_params_has_A, default_meson_mass_guess, meson_mass_equation
using .MottTransition: mott_threshold_mass, mott_gap, is_mott_point

@testset "MesonMass/MottTransition 基础测试" begin
    quark_params = (m=(u=0.3, d=0.31, s=0.5), μ=(u=0.0, d=0.0, s=0.0))
    thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)

    qp = ensure_quark_params_has_A(quark_params, thermo_params)
    G_u = calculate_G_from_A(qp.A.u, qp.m.u)
    G_s = calculate_G_from_A(qp.A.s, qp.m.s)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    @test default_meson_mass_guess(:pi, qp) ≈ qp.m.u + qp.m.d
    @test default_meson_mass_guess(:K, qp) ≈ qp.m.u + qp.m.s

    f = meson_mass_equation(:pi, 0.5, 0.0, 0.0, qp, thermo_params, K_coeffs)
    @test isfinite(real(f)) && isfinite(imag(f))

    threshold = mott_threshold_mass(:K, qp)
    @test threshold ≈ qp.m.u + qp.m.s
    gap = mott_gap(:K, threshold, qp)
    @test is_mott_point(:K, threshold, qp; atol=1e-12)
end