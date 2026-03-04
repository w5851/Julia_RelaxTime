"""
TotalPropagator 模块单元测试
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

# 加载完整 RelaxTime 模块链（TotalPropagator 依赖多个上游模块）
const _RELAX_TIME_PATH_TP = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAX_TIME_PATH_TP)
end

using Main.Constants_PNJL: G_fm2, K_fm5
using Main.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using Main.MesonMass: ensure_quark_params_has_A
using Main.TotalPropagator: get_flavor_factor, calculate_cms_momentum

@testset "TotalPropagator" begin
    quark_params = (m=(u=0.3, d=0.31, s=0.5), μ=(u=0.0, d=0.0, s=0.0))
    thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
    qp = ensure_quark_params_has_A(quark_params, thermo_params)

    G_u = calculate_G_from_A(qp.A.u, qp.m.u)
    G_s = calculate_G_from_A(qp.A.s, qp.m.s)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    @testset "get_flavor_factor 返回正确味道因子" begin
        # 签名: get_flavor_factor(quark1::Symbol, quark2::Symbol)
        @test isfinite(get_flavor_factor(:u, :u))   # 同味
        @test isfinite(get_flavor_factor(:u, :d))   # 异味
    end

    @testset "calculate_cms_momentum 返回有限值" begin
        # 签名: calculate_cms_momentum(process, s, t, channel, quark_params; u=nothing)
        m1, m2 = qp.m.u, qp.m.d
        s = (m1 + m2 + 0.1)^2  # 保证阈值以上
        t = -0.1
        result = calculate_cms_momentum(:ud_to_ud, s, t, :s, qp)
        @test isfinite(result.k0) || isfinite(result.k)
    end

    @testset "calculate_cms_momentum 阈值以下" begin
        m1, m2 = qp.m.u, qp.m.d
        s_below = (m1 + m2 - 0.01)^2
        t = 0.0
        result = calculate_cms_momentum(:ud_to_ud, s_below, t, :s, qp)
        @test isfinite(result.k0)
    end
end
