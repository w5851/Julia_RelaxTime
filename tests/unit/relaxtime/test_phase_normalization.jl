"""Pure algebra and synthetic-profile tests for phase normalization."""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test
using LinearAlgebra: Diagonal

const _RELAXTIME_PATH_PHASE_NORMALIZATION = normpath(joinpath(
    @__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl",
))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAXTIME_PATH_PHASE_NORMALIZATION)
end

using Main.RelaxTime.PhaseNormalization: SCATTERING_PHASE,
                                           SMATRIX_ARGUMENT,
                                           phase_variable,
                                           phase_to_s_matrix,
                                           s_matrix_argument,
                                           s_matrix_to_phase,
                                           phase_measure_factor,
                                           phase_measure,
                                           s_matrix_log_derivative,
                                           s_matrix_density_of_states,
                                           phase_to_s_matrix_profile,
                                           s_matrix_to_phase_profile
using Main.RelaxTime.BUPhaseGates: anchor_phase_high_energy

@testset "Phase variable and measure contract" begin
    @test phase_variable(:delta) === SCATTERING_PHASE
    @test phase_variable(:arg_s) === SMATRIX_ARGUMENT
    @test phase_measure_factor(:delta) ≈ inv(pi)
    @test phase_measure_factor(:s_matrix_argument) ≈ inv(2pi)
    @test phase_measure(:delta) === :d_delta_over_pi
    @test phase_measure(:s_matrix_argument) === :d_arg_s_over_2pi
    @test_throws ArgumentError phase_variable(:propagator)
end

@testset "Scalar S-matrix algebra closes the factor of two" begin
    delta = 0.73
    ddelta = -0.41
    S = phase_to_s_matrix(delta)
    dS = 2im * ddelta * S

    @test abs(S) ≈ 1.0 atol=1e-14
    @test s_matrix_argument(S) ≈ 2delta atol=1e-14
    @test s_matrix_to_phase(S) ≈ delta atol=1e-14
    @test s_matrix_log_derivative(S, dS) ≈ 2ddelta atol=1e-14
    @test s_matrix_density_of_states(S, dS) ≈ ddelta / pi atol=1e-14

    argument = 2delta
    S_argument = phase_to_s_matrix(argument; variable=:s_matrix_argument)
    @test S_argument ≈ S atol=1e-14
    @test s_matrix_to_phase(S_argument; variable=:s_matrix_argument) ≈ argument atol=1e-14
end

@testset "Synthetic profile preserves continuous phase and DMB measure" begin
    omega = collect(range(0.0, 8.0; length=401))
    # Keep the synthetic phase below pi/2 so the scalar S-matrix principal
    # argument is injective; a bound-state pi jump needs an independent
    # endpoint/count convention and cannot be recovered from S alone.
    delta = @. 0.4pi * (1.0 - tanh((omega - 2.0) / 0.3)) / 2.0
    S = phase_to_s_matrix_profile(delta)
    principal = s_matrix_to_phase_profile(S)
    anchored = anchor_phase_high_energy(omega, principal; tail_points=20)

    @test maximum(abs.(anchored.anchored_phase .- delta)) < 1.0e-10
    @test anchored.anchored_phase[end] ≈ 0.0 atol=1.0e-14
    @test maximum(abs.(abs.(S) .- 1.0)) < 1.0e-14

    ddelta = similar(delta)
    ddelta[1] = (delta[2] - delta[1]) / (omega[2] - omega[1])
    ddelta[end] = (delta[end] - delta[end - 1]) / (omega[end] - omega[end - 1])
    @inbounds for i in 2:(length(delta) - 1)
        ddelta[i] = (delta[i + 1] - delta[i - 1]) / (omega[i + 1] - omega[i - 1])
    end
    dS = ComplexF64[2im * ddelta[i] * S[i] for i in eachindex(S)]
    density_measure = [s_matrix_density_of_states(S[i], dS[i]) for i in eachindex(S)]
    @test maximum(abs.(density_measure .- ddelta ./ pi)) < 1.0e-12
end

@testset "Diagonal S-matrix uses the connected trace" begin
    delta = [0.2, -0.37]
    ddelta = [0.11, 0.29]
    S = Diagonal(phase_to_s_matrix.(delta))
    dS = Diagonal(ComplexF64[2im * ddelta[i] * S[i, i] for i in eachindex(delta)])
    @test s_matrix_log_derivative(S, dS) ≈ 2sum(ddelta) atol=1e-14
    @test s_matrix_density_of_states(S, dS) ≈ sum(ddelta) / pi atol=1e-14
end
