"""
AFieldBuilder 模块单元测试
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

const _CONSTANTS_PNJL_PATH_AF = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "constants", "Constants_PNJL.jl"))
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PNJL_PATH_AF)
end
const _GL_PATH_AF = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "integration", "GaussLegendre.jl"))
if !isdefined(Main, :GaussLegendre)
    Base.include(Main, _GL_PATH_AF)
end
const _QUARK_DISTRIBUTION_PATH_AF = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "models", "pnjl_physics", "QuarkDistribution.jl"))
if !isdefined(Main, :PNJLQuarkDistributions)
    Base.include(Main, _QUARK_DISTRIBUTION_PATH_AF)
end
const _QUARK_DISTRIBUTION_ANISO_PATH_AF = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "QuarkDistribution_Aniso.jl"))
if !isdefined(Main, :PNJLQuarkDistributions_Aniso)
    Base.include(Main, _QUARK_DISTRIBUTION_ANISO_PATH_AF)
end
const _OLI_PATH_AF = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "OneLoopIntegrals.jl"))
if !isdefined(Main, :OneLoopIntegrals)
    Base.include(Main, _OLI_PATH_AF)
end
const _OLI_ANISO_PATH_AF = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "OneLoopIntegralsAniso.jl"))
if !isdefined(Main, :OneLoopIntegralsCorrection)
    Base.include(Main, _OLI_ANISO_PATH_AF)
end
const _AFIELD_BUILDER_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "AFieldBuilder.jl"))
if !isdefined(Main, :AFieldBuilder)
    Base.include(Main, _AFIELD_BUILDER_PATH)
end

using Main.AFieldBuilder: build_A_triplet, ensure_quark_params_has_A

@testset "AFieldBuilder" begin
    Main.AFieldBuilder._AUTO_A_WARNED_ANISO[] = false
    Main.AFieldBuilder._AUTO_A_WARNED_ISO[] = false
    quark_params = (m=(u=0.3, d=0.31, s=0.5), μ=(u=0.0, d=0.0, s=0.0))

    @testset "build_A_triplet 各向同性 (ξ=0)" begin
        thermo_iso = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
        A_tri = build_A_triplet(quark_params, thermo_iso)
        @test haskey(A_tri, :u) && haskey(A_tri, :d) && haskey(A_tri, :s)
        @test all(isfinite, (A_tri.u, A_tri.d, A_tri.s))
    end

    @testset "build_A_triplet 各向异性 (ξ>0)" begin
        thermo_aniso = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.3)
        A_tri = build_A_triplet(quark_params, thermo_aniso)
        @test all(isfinite, (A_tri.u, A_tri.d, A_tri.s))
    end

    @testset "ensure_quark_params_has_A 自动补充 A 字段" begin
        thermo = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
        qp = ensure_quark_params_has_A(quark_params, thermo)
        @test haskey(qp, :A)
        @test haskey(qp.A, :u) && haskey(qp.A, :d) && haskey(qp.A, :s)
    end

    @testset "已有 A 时不重复计算" begin
        thermo = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
        qp_with_A = (m=quark_params.m, μ=quark_params.μ, A=(u=1.0, d=1.0, s=1.0))
        qp2 = ensure_quark_params_has_A(qp_with_A, thermo)
        @test qp2.A.u ≈ 1.0
        @test qp2.A.d ≈ 1.0
        @test qp2.A.s ≈ 1.0
    end

    @testset "ξ≠0 自动补 A 仅告警一次" begin
        thermo = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.3)
        Main.AFieldBuilder._AUTO_A_WARNED_ANISO[] = false
        @test_logs (:warn, r"auto-building anisotropic A via A_aniso") ensure_quark_params_has_A(quark_params, thermo)
        @test_logs ensure_quark_params_has_A(quark_params, thermo)
    end

    @testset "可通过环境变量关闭自动补 A 告警" begin
        thermo = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.3)
        withenv("RELAXTIME_AUTO_A_WARN" => "0") do
            Main.AFieldBuilder._AUTO_A_WARNED_ANISO[] = false
            @test_logs ensure_quark_params_has_A(quark_params, thermo)
        end
    end
end
