using Test

include(joinpath(@__DIR__, "..", "..", "..", "scripts", "analysis", "pnjl_bifurcation_spike", "matched_pilot_support.jl"))
using .MatchedPALCPilotSupport: target_seeds, NewtonWork, observe_newton!, finish_attempt!

@testset "PALC target correction seeds preserve crossings" begin
    folded = [(rho=0.0, solution=[0.0]), (rho=2.0, solution=[2.0]),
        (rho=0.0, solution=[4.0])]
    @test target_seeds(folded, 1.0) == [[1.0], [3.0]]
    @test target_seeds(folded, 2.0) == [[2.0]]
    @test target_seeds(folded, 3.0) == [[2.0]]
    @test isempty(target_seeds(NamedTuple[], 1.0))
end

@testset "PALC Newton telemetry includes failed retries once" begin
    work = NewtonWork()
    history = [1.0]
    state(step, residual, h=history) = (; step, residual, residuals=h, options=(tol=1e-8,))
    observe_newton!(work, state(0, 1.0))
    observe_newton!(work, state(1, 0.1))
    observe_newton!(work, state(1, 0.1)) # final callback repeats the last step
    next_history = [0.1]
    observe_newton!(work, state(0, 0.1, next_history))
    observe_newton!(work, state(1, 1e-10, next_history))
    observe_newton!(work, state(1, 1e-10, next_history))
    finish_attempt!(work)
    finish_attempt!(work)
    @test work.attempts == 2
    @test work.iterations == 2
    @test work.nonconverged == 1
end
