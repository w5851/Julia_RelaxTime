using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
const RUNNER = joinpath(PROJECT_ROOT, "scripts", "analysis",
    "pnjl_maxwell_endpoint_refinement.jl")
if !isdefined(Main, :MaxwellEndpointRefinementLoaded)
    include(RUNNER)
    const MaxwellEndpointRefinementLoaded = true
end
const REPLAY = Main.MaxwellEndpointCandidateFeasibility

function _small_curve(rho; mu=[2.0, 4.0, 2.0, 1.0, 3.0, 5.0])
    points = [REPLAY.CurvePoint(-0.5, 5.0, "synthetic", Float64(r), Float64(m),
        0.0, 0, "base_grid") for (r, m) in zip(rho, mu)]
    REPLAY.CurveData(points, Float64.(rho), Float64.(mu))
end

@testset "Maxwell endpoint refinement routing" begin
    curve = _small_curve([0.0, 0.001, 0.002, 0.003, 0.004, 0.005])
    result = REPLAY.strict_candidate(curve)
    @test result.status == :first_order
    bracket = Main._bracket(curve, result)
    @test bracket !== nothing
    @test !bracket.positive
    @test bracket.left == 0.0
    target = Main._target_for_next(curve, result, Set{Float64}())
    @test target == 0.0005

    targeted_curve = _small_curve([0.0001, 0.0005, 0.001, 0.002, 0.003, 0.004])
    targeted_result = (status=:first_order,
        candidate=(crossings=[0.0002, 0.0015, 0.0035],))
    targeted_bracket = Main._bracket(targeted_curve, targeted_result)
    @test targeted_bracket !== nothing
    @test targeted_bracket.positive
    @test targeted_bracket.width <= 0.001
end
