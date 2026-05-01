"""
MesonDensity 模块单元测试（稳定粒子极限）。
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

const _RELAXTIME_PATH_MD = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAXTIME_PATH_MD)
end

using Main.RelaxTime.MesonDensity: DEFAULT_MESON_DENSITY_Q_NODES,
                                   bose_distribution,
                                   meson_degeneracy,
                                   stable_meson_number_density,
                                   stable_kpi_ratio,
                                   stable_kpi_scan

const ZETA3 = 1.2020569031595942

@testset "MesonDensity 简并因子与玻色分布" begin
    @test meson_degeneracy(:pi) == 3
    @test meson_degeneracy(:K) == 4
    @test meson_degeneracy(:pi; charge_resolved=true) == 1
    @test meson_degeneracy(:K; charge_resolved=true) == 1

    @test bose_distribution(1.0, 0.0, 0.2) > 0.0
    @test bose_distribution(1.0, 0.0, 0.0) == 0.0
    @test_throws ArgumentError bose_distribution(0.5, 0.5, 0.2)
end

@testset "MesonDensity 稳定粒子极限" begin
    T = 0.25
    n_pi = stable_meson_number_density(0.14, T; degeneracy=3, num_q_nodes=DEFAULT_MESON_DENSITY_Q_NODES)
    n_K = stable_meson_number_density(0.49, T; degeneracy=4, num_q_nodes=DEFAULT_MESON_DENSITY_Q_NODES)

    @test n_pi > 0.0
    @test n_K > 0.0
    @test n_pi > n_K

    n_low = stable_meson_number_density(0.14, 0.15; degeneracy=3)
    n_high = stable_meson_number_density(0.14, 0.20; degeneracy=3)
    @test n_high > n_low

    # 质量为零、μ=0 时的解析结果：n = d * ζ(3) / π² * T³
    n_massless = stable_meson_number_density(0.0, 0.3; degeneracy=2, μ=-1e-6, qmax=12.0, num_q_nodes=512)
    n_expected = 2.0 * ZETA3 / π^2 * 0.3^3
    @test n_massless ≈ n_expected rtol=2e-3

    @test_throws ArgumentError stable_meson_number_density(0.14, 0.2; μ=0.14)
end

@testset "MesonDensity K/π 扫描" begin
    ratio = stable_kpi_ratio(0.14, 0.49, 0.20)
    @test 0.0 < ratio < 1.0

    scan = stable_kpi_scan([0.16, 0.18, 0.20]; m_pi=0.14, m_K=0.49)
    @test length(scan.temperatures) == 3
    @test length(scan.n_pi) == 3
    @test length(scan.n_K) == 3
    @test length(scan.kpi_ratio) == 3
    @test all(diff(scan.n_pi) .> 0.0)
    @test all(diff(scan.n_K) .> 0.0)
end
