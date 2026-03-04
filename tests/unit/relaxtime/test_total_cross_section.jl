"""
TotalCrossSection 模块单元测试
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

# 加载完整 RelaxTime 模块链
const _RELAX_TIME_PATH_TCS = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAX_TIME_PATH_TCS)
end

using Main.Constants_PNJL: G_fm2, K_fm5
using Main.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using Main.MesonMass: ensure_quark_params_has_A
using Main.TotalCrossSection: calculate_t_bounds, calculate_final_state_energies, final_state_blocking_factor

@testset "TotalCrossSection" begin
    quark_params = (m=(u=0.3, d=0.31, s=0.5), μ=(u=0.0, d=0.0, s=0.0))
    thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
    qp = ensure_quark_params_has_A(quark_params, thermo_params)

    m1 = qp.m.u
    m2 = qp.m.d
    m3 = qp.m.u  # 散射后味道
    m4 = qp.m.d
    s = (m1 + m2 + 0.2)^2

    @testset "calculate_t_bounds 返回合理范围" begin
        t_min, t_max = calculate_t_bounds(s, m1, m2, m3, m4)
        @test isfinite(t_min)
        @test isfinite(t_max)
        @test t_min ≤ t_max
    end

    @testset "calculate_final_state_energies 返回有限能量" begin
        E3, E4 = calculate_final_state_energies(s, m3, m4)
        @test E3 > 0
        @test E4 > 0
        @test isfinite(E3)
        @test isfinite(E4)
    end

    @testset "final_state_blocking_factor 在 [0,1] 范围" begin
        # 签名: (flavor::Symbol, E, m, μ, T, Φ, Φbar, ξ, cosθ)
        f = final_state_blocking_factor(:u, 0.5, m3, 0.0, 0.15, 0.5, 0.5, 0.0, 0.5)
        @test 0.0 ≤ f ≤ 1.0
    end
end
