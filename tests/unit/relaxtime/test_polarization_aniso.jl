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
const _OLI_PATH_PA = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "OneLoopIntegrals.jl"))
if !isdefined(Main, :OneLoopIntegrals)
    Base.include(Main, _OLI_PATH_PA)
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
end
