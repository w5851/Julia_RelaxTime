"""
MottTransition 模块单元测试（轻量验证）。
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

const _CONSTANTS_PNJL_PATH_MT = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "constants", "Constants_PNJL.jl"))
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PNJL_PATH_MT)
end
const _EFFECTIVE_COUPLINGS_PATH_MT = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "EffectiveCouplings.jl"))
if !isdefined(Main, :EffectiveCouplings)
    Base.include(Main, _EFFECTIVE_COUPLINGS_PATH_MT)
end
const _MESON_MASS_PATH_MT = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "MesonMass.jl"))
if !isdefined(Main, :MesonMass)
    Base.include(Main, _MESON_MASS_PATH_MT)
end
const _MOTT_TRANSITION_PATH_MT = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "MottTransition.jl"))
if !isdefined(Main, :MottTransition)
    Base.include(Main, _MOTT_TRANSITION_PATH_MT)
end

using Main.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using Main.Constants_PNJL: G_fm2, K_fm5
using Main.MesonMass: ensure_quark_params_has_A
using Main.MottTransition: mott_threshold_mass, mott_gap, is_mott_point, mott_threshold_masses, mott_gaps

@testset "MottTransition 基础测试" begin
    quark_params = (m=(u=0.3, d=0.31, s=0.5), μ=(u=0.0, d=0.0, s=0.0))
    thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)

    qp = ensure_quark_params_has_A(quark_params, thermo_params)

    threshold = mott_threshold_mass(:K, qp)
    @test threshold ≈ qp.m.u + qp.m.s
    gap = mott_gap(:K, threshold, qp)
    @test is_mott_point(:K, threshold, qp; atol=1e-12)

    # Mixed-meson threshold API
    thr_eta = mott_threshold_masses(:eta, qp)
    @test thr_eta.uu ≈ 2.0 * qp.m.u
    @test thr_eta.ss ≈ 2.0 * qp.m.s
    @test thr_eta.min ≈ min(thr_eta.uu, thr_eta.ss)

    gaps_eta = mott_gaps(:eta, thr_eta.min, qp)
    @test isfinite(gaps_eta.uu) && isfinite(gaps_eta.ss) && isfinite(gaps_eta.min)
end
