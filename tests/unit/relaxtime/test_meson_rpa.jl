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

    @testset "charged ladder vertices derive scalar/matrix normalization" begin
        lambda1 = @SMatrix ComplexF64[
            0 1 0
            1 0 0
            0 0 0
        ]
        lambda2 = @SMatrix ComplexF64[
            0 -im 0
            im 0 0
            0 0 0
        ]
        t_plus = (lambda1 + im * lambda2) / sqrt(2)
        t_minus = (lambda1 - im * lambda2) / sqrt(2)

        bubble_ud = 0.12 + 0.03im
        bubble_du = 0.08 - 0.02im
        ordered_bubbles = @SMatrix ComplexF64[
            0 bubble_ud 0
            bubble_du 0 0
            0 0 0
        ]
        project(left, right) = sum(
            left[i, j] * right[j, i] * ordered_bubbles[i, j]
            for i in 1:3, j in 1:3
        )

        # T_+ = sqrt(2) E_ud and T_- = sqrt(2) E_du.  The factor two
        # therefore comes from the two charged flavor vertices, not from an
        # assumed conversion between two already-computed polarizations.
        @test t_plus[1, 2] ≈ sqrt(2)
        @test t_minus[2, 1] ≈ sqrt(2)
        @test count(!iszero, t_plus) == 1
        @test count(!iszero, t_minus) == 1
        @test project(t_plus, t_minus) ≈ 2 * bubble_ud
        @test project(t_minus, t_plus) ≈ 2 * bubble_du
        @test project(lambda1, lambda1) ≈ bubble_ud + bubble_du
        @test project(lambda1, lambda2) ≈ im * (bubble_ud - bubble_du)

        coupling = 0.37
        charged = 2 * coupling / (1 - 4 * coupling * bubble_ud)
        matrix_charge_basis = 2 * coupling / (1 - 2 * coupling * project(t_plus, t_minus))
        unconverted_matrix = 2 * coupling / (1 - 2 * coupling * bubble_ud)

        @test charged ≈ matrix_charge_basis
        @test !isapprox(charged, unconverted_matrix; rtol=1e-12, atol=1e-14)
    end

    @testset "chiral-limit Goldstone identity fixes the factor four" begin
        Nc = 3.0
        A_u = -2.4
        A_s = -1.7
        m_u = 1.8
        m_s = 2.6
        K_six = 0.41

        H_s = -Nc * m_s * A_s / (4 * pi^2)
        K12 = -pi^2 / (Nc * A_u)
        G = K12 - K_six * H_s / 2

        gap_rhs = -G * Nc * m_u * A_u / pi^2 +
                  K_six * Nc^2 * m_u * A_u * m_s * A_s / (8 * pi^4)
        pi_rehberg = -Nc * A_u / (4 * pi^2)
        charged_inverse = 1 - 4 * K12 * pi_rehberg
        matrix_inverse = 1 - 2 * K12 * (2 * pi_rehberg)
        wrong_inverse = 1 - 2 * K12 * pi_rehberg

        @test gap_rhs ≈ m_u
        @test charged_inverse ≈ 0.0 atol=1e-14
        @test matrix_inverse ≈ 0.0 atol=1e-14
        @test !isapprox(wrong_inverse, 0.0; atol=1e-14)
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
