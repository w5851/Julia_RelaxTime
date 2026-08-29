"""Unit tests for the Phase-1 full KMT interaction-kernel backend."""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test
using LinearAlgebra: I
using StaticArrays

const _ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _RELAXTIME_PATH = joinpath(_ROOT, "src", "relaxtime", "RelaxTime.jl")
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAXTIME_PATH)
end

using Main.RelaxTime.MesonInteractionKernel
using Main.RelaxTime.EffectiveCouplings: calculate_effective_couplings

@testset "MesonInteractionKernel Phase 1" verbose=true begin
    G = 1.25
    K = 0.40

    @testset "结构与输入适配" begin
        kernel = build_full_kmt_interaction(SVector(-0.30, -0.20, -0.10); G=G, K=K)
        @test kernel isa FullKMTInteraction{Float64}
        @test Main.MesonInteractionKernel === Main.RelaxTime.MesonInteractionKernel
        @test kernel.neutral_basis == (:lambda0, :lambda3, :lambda8)
        @test size(neutral_coupling_matrix(kernel, :P)) == (3, 3)
        @test charged_couplings(kernel, :S) isa ChargedKMTChannels{Float64}

        @test build_full_kmt_interaction([-0.30, -0.20, -0.10]; G=G, K=K).phi == kernel.phi
        @test build_full_kmt_interaction((-0.30, -0.20, -0.10); G=G, K=K).phi == kernel.phi
        @test build_full_kmt_interaction((u=-0.30, d=-0.20, s=-0.10); G=G, K=K).phi == kernel.phi

        integer_kernel = build_full_kmt_interaction((1, 2, 3); G=1, K=1)
        @test integer_kernel isa FullKMTInteraction{Float64}
        @test all(isfinite, integer_kernel.neutral_P)

        equilibrium = (converged=true, x_state=SVector(-0.30, -0.20, -0.10, 0.4, 0.4))
        adapted = build_full_kmt_interaction_from_equilibrium(equilibrium; G=G, K=K)
        @test adapted.phi == kernel.phi
        @test adapted.neutral_P == kernel.neutral_P
    end

    @testset "不对称凝聚产生 0/3/8 交叉项" begin
        u, d, strange = -0.30, -0.20, -0.10
        kernel = build_full_kmt_interaction((u, d, strange); G=G, K=K)
        expected_K03 = K / (2 * sqrt(6)) * (u - d)
        expected_K08 = -sqrt(2) * K / 12 * (u + d - 2 * strange)
        expected_K38 = -K / (2 * sqrt(3)) * (u - d)

        @test get_Kab(kernel, 0, 3, :P) ≈ expected_K03
        @test get_Kab(kernel, 3, 0, :P) ≈ expected_K03
        @test get_Kab(kernel, 0, 8, :P) ≈ expected_K08
        @test get_Kab(kernel, 3, 8, :P) ≈ expected_K38
        @test get_Kab(kernel, 8, 3, :P) ≈ expected_K38
        @test charged_coupling(kernel, :K45, :P) ≈ G - K * d / 2
        @test charged_coupling(kernel, :K67, :P) ≈ G - K * u / 2
        @test charged_coupling(kernel, :K45, :P) != charged_coupling(kernel, :K67, :P)

        @test get_Kab(kernel, 0, 3, :S) ≈ -expected_K03
        @test get_Kab(kernel, 0, 8, :S) ≈ -expected_K08
        @test get_Kab(kernel, 3, 8, :S) ≈ -expected_K38
    end

    @testset "同位旋对称极限与旧 0/8 API 兼容" begin
        phi_l, phi_s = -0.30, -0.12
        kernel = build_full_kmt_interaction((phi_l, phi_l, phi_s); G=G, K=K)
        old = calculate_effective_couplings(G, K, -phi_l, -phi_s)
        P, S = kernel.neutral_P, kernel.neutral_S

        @test P[1, 2] == 0.0
        @test P[2, 3] == 0.0
        @test charged_coupling(kernel, :K45, :P) == charged_coupling(kernel, :K67, :P)
        @test charged_coupling(kernel, :K45, :S) == charged_coupling(kernel, :K67, :S)

        @test P[1, 1] ≈ old.K0_plus
        @test P[2, 2] ≈ old.K123_plus
        @test P[3, 3] ≈ old.K8_plus
        @test P[1, 3] ≈ old.K08_plus
        @test charged_coupling(kernel, :K45, :P) ≈ old.K4567_plus

        @test S[1, 1] ≈ old.K0_minus
        @test S[2, 2] ≈ old.K123_minus
        @test S[3, 3] ≈ old.K8_minus
        @test S[1, 3] ≈ old.K08_minus
        @test charged_coupling(kernel, :K45, :S) ≈ old.K4567_minus
    end

    @testset "KMT 关闭与 P/S 互补关系" begin
        kernel = build_full_kmt_interaction((-0.30, -0.20, -0.10); G=G, K=0.0)
        @test kernel.neutral_P == G .* I(3)
        @test kernel.neutral_S == G .* I(3)
        @test kernel.charged_P == ChargedKMTChannels(G, G, G)
        @test kernel.charged_S == ChargedKMTChannels(G, G, G)

        kernel = build_full_kmt_interaction((-0.30, -0.20, -0.10); G=G, K=K)
        @test kernel.neutral_P + kernel.neutral_S ≈ 2G .* I(3)
        @test kernel.charged_P.K12 + kernel.charged_S.K12 ≈ 2G
        @test kernel.charged_P.K45 + kernel.charged_S.K45 ≈ 2G
        @test kernel.charged_P.K67 + kernel.charged_S.K67 ≈ 2G
    end

    @testset "参数和查询校验" begin
        @test_throws ArgumentError build_full_kmt_interaction([1.0, 2.0]; G=G, K=K)
        @test_throws ArgumentError build_full_kmt_interaction(42; G=G, K=K)
        @test_throws ArgumentError build_full_kmt_interaction((NaN, 0.0, 0.0); G=G, K=K)
        @test_throws ArgumentError build_full_kmt_interaction((0.0, 0.0, 0.0); G=Inf, K=K)
        kernel = build_full_kmt_interaction((0.0, 0.0, 0.0); G=G, K=K)
        @test_throws ArgumentError neutral_coupling_matrix(kernel, :unknown)
        @test_throws ArgumentError charged_coupling(kernel, :K99, :P)
        @test_throws ArgumentError get_Kab(kernel, 1, 3, :P)
        @test_throws ArgumentError build_full_kmt_interaction_from_equilibrium((state=[1, 2, 3],); G=G, K=K)
        @test_throws ArgumentError build_full_kmt_interaction_from_equilibrium((x_state=[1, 2],); G=G, K=K)
    end
end
