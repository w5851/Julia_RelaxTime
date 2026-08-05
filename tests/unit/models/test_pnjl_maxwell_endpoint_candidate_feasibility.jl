using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
const REPLAY_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "analysis",
    "pnjl_maxwell_endpoint_candidate_feasibility.jl")
if !isdefined(Main, :MaxwellEndpointCandidateFeasibility)
    include(REPLAY_SCRIPT)
end
const REPLAY = Main.MaxwellEndpointCandidateFeasibility

function _curve(mu; rho=collect(0.0:(length(mu) - 1)))
    points = [REPLAY.CurvePoint(-0.5, 5.0, "synthetic", Float64(r), Float64(m),
        0.0, 0, "grid") for (r, m) in zip(rho, mu)]
    REPLAY.CurveData(points, Float64.(rho), Float64.(mu))
end

@testset "Maxwell endpoint/candidate feasibility replay" begin
    @testset "two-crossing endpoint jump is invalid" begin
        curve = _curve([2.0, 4.0, 2.0, 1.0, 3.0, 5.0])
        sres = Models.detect_s_shape(curve.mu, curve.rho)
        low = REPLAY._area_at(1.5, curve, sres)
        high = REPLAY._area_at(2.5, curve, sres)
        @test !low.valid
        @test length(low.crossings) == 2
        @test high.valid
        @test length(high.crossings) == 3
    end

    @testset "unique three-crossing candidate is found" begin
        curve = _curve([2.0, 4.0, 2.0, 1.0, 3.0, 5.0];
            rho=[0.0, 0.001, 0.002, 0.003, 0.004, 0.005])
        result = REPLAY.strict_candidate(curve)
        @test result.status == :first_order
        @test result.reason == "unique_three_crossing_candidate"
        @test length(result.roots) == 1
        @test isapprox(result.candidate.mu, 2.357143; atol=2e-4)
        @test result.endpoint_dependent
    end

    @testset "crossing vertices are deduplicated" begin
        curve = _curve([0.0, 2.0, 1.0, 0.0, 2.0, 3.0])
        crossings = REPLAY._crossings(2.0, curve)
        @test crossings == [1.0, 4.0]
    end

    @testset "max iteration is a failure, not a converged result" begin
        curve = _curve([2.0, 4.0, 2.0, 1.0, 3.0, 5.0])
        sres = Models.detect_s_shape(curve.mu, curve.rho)
        probe_a = REPLAY._area_at(2.35, curve, sres)
        probe_b = REPLAY._area_at(2.36, curve, sres)
        solved = REPLAY._bisect(curve, sres, 2.35, 2.36,
            probe_a.area, probe_b.area; max_iter=1)
        @test solved !== nothing
        @test !solved.converged
        @test solved.reason == "solver_tolerance_not_met"
    end

    @testset "monotone curve is not a Maxwell candidate" begin
        curve = _curve(collect(0.0:1.0:15.0))
        result = REPLAY.strict_candidate(curve)
        @test result.status == :monotone
        @test isempty(result.roots)
    end

    @testset "endpoint policy never creates duplicate rho zero points" begin
        curve = _curve([2.0, 4.0, 2.0, 1.0, 3.0, 5.0];
            rho=[0.0, 0.001, 0.002, 0.003, 0.004, 0.005])
        @test count(==(0.0), curve.rho) == 1
        result = REPLAY.strict_candidate(curve)
        @test result.endpoint_dependent
    end
end
