"""
MesonMass 模块单元测试（轻量验证）。
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

const _CONSTANTS_PNJL_PATH_MM = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "constants", "Constants_PNJL.jl"))
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PNJL_PATH_MM)
end
const _EFFECTIVE_COUPLINGS_PATH_MM = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "EffectiveCouplings.jl"))
if !isdefined(Main, :EffectiveCouplings)
    Base.include(Main, _EFFECTIVE_COUPLINGS_PATH_MM)
end
const _MESON_MASS_PATH_MM = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "MesonMass.jl"))
if !isdefined(Main, :MesonMass)
    Base.include(Main, _MESON_MASS_PATH_MM)
end

using Main.Constants_PNJL: G_fm2, K_fm5
using Main.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using Main.MesonMass: ensure_quark_params_has_A, default_meson_mass_guess, meson_mass_equation

@testset "MesonMass 基础测试" begin
    quark_params = (m=(u=0.3, d=0.31, s=0.5), μ=(u=0.0, d=0.0, s=0.0))
    thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)

    qp = ensure_quark_params_has_A(quark_params, thermo_params)
    G_u = calculate_G_from_A(qp.A.u, qp.m.u)
    G_s = calculate_G_from_A(qp.A.s, qp.m.s)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    @test default_meson_mass_guess(:pi, qp) ≈ qp.m.u + qp.m.d
    @test default_meson_mass_guess(:K, qp) ≈ qp.m.u + qp.m.s
    @test default_meson_mass_guess(:eta, qp) ≈ 2.0 * qp.m.u
    @test default_meson_mass_guess(:eta_prime, qp) ≈ 2.0 * qp.m.s

    f = meson_mass_equation(:pi, 0.5, 0.0, 0.0, qp, thermo_params, K_coeffs)
    @test isfinite(real(f)) && isfinite(imag(f))
end

@testset "MesonMass 混合通道一致性" begin
    # Construct an (approx) SU(3)-symmetric point so Π_uu == Π_ss and K08 == 0
    quark_params = (m=(u=0.35, d=0.35, s=0.35), μ=(u=0.0, d=0.0, s=0.0))
    thermo_params = (T=0.12, Φ=0.6, Φbar=0.6, ξ=0.0)

    qp = ensure_quark_params_has_A(quark_params, thermo_params)
    G_u = calculate_G_from_A(qp.A.u, qp.m.u)
    G_s = calculate_G_from_A(qp.A.s, qp.m.s)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    # Access internal helper for matrix elements
    elems = MesonMass._mixed_matrix_elements(:P, 0.7, 0.0, 0.0, qp, thermo_params, K_coeffs)
    @test abs(elems.M08) ≤ 1e-10

    f_eta = meson_mass_equation(:eta, 0.7, 0.0, 0.0, qp, thermo_params, K_coeffs)
    f_etap = meson_mass_equation(:eta_prime, 0.7, 0.0, 0.0, qp, thermo_params, K_coeffs)
    detM = elems.M00 * elems.M88 - elems.M08^2
    @test f_eta * f_etap ≈ 4.0 * detM atol=1e-10 rtol=1e-10
end
