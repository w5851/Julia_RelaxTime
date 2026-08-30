"""
PolarizationAniso 模块单元测试
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

const _CONSTANTS_PNJL_PATH_PA = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "constants", "Constants_PNJL.jl"))
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PNJL_PATH_PA)
end
const _GAUSS_LEGENDRE_PATH_PA = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "integration", "GaussLegendre.jl"))
if !isdefined(Main, :GaussLegendre)
    Base.include(Main, _GAUSS_LEGENDRE_PATH_PA)
end
const _QUARK_DISTRIBUTION_PATH_PA = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "models", "pnjl_physics", "QuarkDistribution.jl"))
if !isdefined(Main, :PNJLQuarkDistributions)
    Base.include(Main, _QUARK_DISTRIBUTION_PATH_PA)
end
const _OLI_PATH_PA = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "OneLoopIntegrals.jl"))
if !isdefined(Main, :OneLoopIntegrals)
    Base.include(Main, _OLI_PATH_PA)
end
const _QUARK_DISTRIBUTION_ANISO_PATH_PA = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "QuarkDistribution_Aniso.jl"))
if !isdefined(Main, :PNJLQuarkDistributions_Aniso)
    Base.include(Main, _QUARK_DISTRIBUTION_ANISO_PATH_PA)
end
const _OLI_ANISO_PATH_PA = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "OneLoopIntegralsAniso.jl"))
if !isdefined(Main, :OneLoopIntegralsCorrection)
    Base.include(Main, _OLI_ANISO_PATH_PA)
end
const _POLARIZATION_ANISO_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "PolarizationAniso.jl"))
if !isdefined(Main, :PolarizationAniso)
    Base.include(Main, _POLARIZATION_ANISO_PATH)
end

using Main.PolarizationAniso: polarization_aniso, polarization_with_width, polarization_complex

@testset "PolarizationAniso" begin
    # 基本参数
    m1, m2 = 0.3, 0.31
    μ1, μ2 = 0.0, 0.0
    T, Φ, Φbar = 0.15, 0.5, 0.5
    ξ = 0.3  # 各向异性参数
    A1, A2 = 1.0, 1.0
    num_s = 24

    @testset "polarization_aniso 返回有限值" begin
        for ch in (:P, :S)
            re, im = polarization_aniso(ch, 0.5, 0.0, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, num_s)
            @test isfinite(re)
            @test isfinite(im)
        end
    end

    @testset "ξ=0 时退化为各向同性" begin
        re_aniso, im_aniso = polarization_aniso(:P, 0.5, 0.0, m1, m2, μ1, μ2, T, Φ, Φbar, 0.0, A1, A2, num_s)
        @test isfinite(re_aniso)
        @test isfinite(im_aniso)
    end

    @testset "polarization_with_width 返回有限值" begin
        re, im = polarization_with_width(:P, 0.5, 0.01, 0.0, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, num_s)
        @test isfinite(re)
        @test isfinite(im)
    end

    @testset "polarization_complex 返回有限值" begin
        p0 = ComplexF64(0.5, 0.01)
        re, im = polarization_complex(:P, p0, 0.0, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, num_s)
        @test isfinite(re)
        @test isfinite(im)
    end

    @testset "finite-mu flavor order and legacy symmetrization stay distinct" begin
        ordered_us = polarization_aniso(
            :P, 0.5, 0.0, 0.30, 0.55, 0.12, -0.04,
            T, Φ, Φbar, 0.0, 1.0, 0.8, 0,
        )
        ordered_su = polarization_aniso(
            :P, 0.5, 0.0, 0.55, 0.30, -0.04, 0.12,
            T, Φ, Φbar, 0.0, 0.8, 1.0, 0,
        )
        legacy_us = polarization_aniso(
            :P, 0.5, 0.0, 0.30, 0.55, 0.12, -0.04,
            T, Φ, Φbar, 0.0, 1.0, 0.8, 1,
        )

        @test all(isfinite, ordered_us)
        @test all(isfinite, ordered_su)
        @test all(isfinite, legacy_us)
        @test any(!isapprox(a, b; rtol=1e-8, atol=1e-10) for (a, b) in zip(ordered_us, ordered_su))
        @test any(!isapprox(a, b; rtol=1e-8, atol=1e-10) for (a, b) in zip(ordered_us, legacy_us))
    end
end
