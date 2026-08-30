"""Unit tests for the source-aligned neutral three-flavor RPA algebra."""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test
using LinearAlgebra: det, I
using StaticArrays

const _ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _RELAXTIME_PATH = joinpath(_ROOT, "src", "relaxtime", "RelaxTime.jl")
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAXTIME_PATH)
end

using Main.RelaxTime.MesonInteractionKernel: build_full_kmt_interaction, neutral_coupling_matrix
using Main.RelaxTime.MesonRPA

@testset "MesonRPA neutral algebra" verbose=true begin
    G = 1.25
    K = 0.40
    kernel = build_full_kmt_interaction((-0.30, -0.20, -0.10); G=G, K=K)

    @testset "flavor bubbles to (0,3,8) matrix" begin
        u, d, s = 0.11 + 0.01im, 0.07 - 0.02im, 0.03 + 0.04im
        pi = neutral_polarization_matrix((u, d, s))

        @test size(pi) == (3, 3)
        @test pi[1, 1] ≈ (2 / 3) * (u + d + s)
        @test pi[2, 2] ≈ u + d
        @test pi[3, 3] ≈ (u + d + 4s) / 3
        @test pi[1, 2] ≈ sqrt(6) * (u - d) / 3
        @test pi[1, 3] ≈ sqrt(2) * (u + d - 2s) / 3
        @test pi[2, 3] ≈ sqrt(3) * (u - d) / 3
        @test pi == transpose(pi)

        @test neutral_polarization_matrix((u=u, d=d, s=s)) == pi
        @test_throws ArgumentError neutral_polarization_matrix((u, d))
        @test_throws ArgumentError neutral_polarization_matrix((u, d, Inf + 0im))
    end

    @testset "source RPA matrix order" begin
        bubbles = SVector(
            0.11 + 0.01im,
            0.07 - 0.02im,
            0.03 + 0.04im,
        )
        pi = neutral_polarization_matrix(bubbles)
        Kp = neutral_coupling_matrix(kernel, :P)
        denominator = neutral_rpa_inverse_matrix(kernel, pi; channel=:P)
        propagator = neutral_rpa_propagator(kernel, pi; channel=:P)

        K_dense = ComplexF64.(Matrix(Kp))
        pi_dense = ComplexF64.(Matrix(pi))
        expected_denominator = Matrix{ComplexF64}(I, 3, 3) - 2 .* K_dense * pi_dense
        expected_propagator = 2 .* K_dense * inv(expected_denominator)

        @test Matrix(denominator) ≈ expected_denominator
        @test Matrix(propagator) ≈ expected_propagator
        @test neutral_rpa_determinant(kernel, pi; channel=:P) ≈ det(expected_denominator)
        @test all(isfinite, real.(propagator))
        @test all(isfinite, imag.(propagator))
    end

    @testset "isospin-symmetric decoupling" begin
        symmetric_kernel = build_full_kmt_interaction((-0.30, -0.30, -0.10); G=G, K=K)
        pi = neutral_polarization_matrix((0.11 + 0.01im, 0.11 + 0.01im, 0.03 + 0.04im))
        propagator = neutral_rpa_propagator(symmetric_kernel, pi; channel=:P)

        @test symmetric_kernel.neutral_P[1, 2] == 0.0
        @test symmetric_kernel.neutral_P[2, 3] == 0.0
        @test pi[1, 2] == 0.0 + 0.0im
        @test pi[2, 3] == 0.0 + 0.0im
        @test propagator[1, 2] ≈ 0.0 + 0.0im
        @test propagator[2, 3] ≈ 0.0 + 0.0im
        @test propagator[2, 1] ≈ 0.0 + 0.0im
        @test propagator[3, 2] ≈ 0.0 + 0.0im
    end

    @testset "KMT-off diagonal polarization remains supported" begin
        free_kernel = build_full_kmt_interaction((-0.30, -0.20, -0.10); G=G, K=0.0)
        pi = neutral_polarization_matrix((0.11 + 0im, 0.07 + 0im, 0.03 + 0im))
        propagator = neutral_rpa_propagator(free_kernel, pi; channel=:P)
        @test free_kernel.neutral_P == G .* I(3)
        @test propagator[1, 2] != 0.0 + 0.0im
        @test propagator[2, 3] != 0.0 + 0.0im
    end

    @testset "channel and matrix validation" begin
        pi = neutral_polarization_matrix((0.1, 0.08, 0.03))
        @test neutral_rpa_propagator(kernel, pi; channel=:S) isa SMatrix{3,3,ComplexF64,9}
        @test_throws ArgumentError neutral_rpa_propagator(kernel, pi; channel=:unknown)
        @test_throws ArgumentError neutral_rpa_propagator(kernel, zeros(2, 2); channel=:P)
        @test_throws ArgumentError neutral_rpa_propagator(kernel, fill(Inf + 0im, 3, 3); channel=:P)
    end
end
