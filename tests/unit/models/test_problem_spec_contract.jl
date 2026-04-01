using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Constants_PNJL)
    include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
end

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "problem spec contract" begin
    @test isdefined(Models, :ProblemSpec)
    @test isdefined(Models, :build_problem_spec)

    @testset "mode dimensions" begin
        modes = [
            Models.FixedMu(),
            Models.FixedRho(1.0),
            Models.FixedAsymmetricRho(1.0, 0.876, 0.0),
            Models.FixedEntropy(0.5),
            Models.FixedSigma(8.0),
        ]

        for mode in modes
            spec = Models.build_problem_spec(mode)
            @test spec isa Models.ProblemSpec
            @test spec.mode == mode
            @test spec.x_dim == Models.state_dim(mode)
            @test spec.theta_dim == Models.param_dim(mode)
            @test spec.residual! isa Function
            @test spec.forward_solve isa Function
            @test spec.conditions isa Function
            @test spec.unpack_solution isa Function
            @test spec.postprocess isa Function
            @test spec.hard_rules isa AbstractVector
            @test spec.selector isa Function
        end
    end
end
