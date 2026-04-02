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

    @testset "fixedrho spec conditions and forward solve are wired" begin
        model = Models.create_model(:PNJL)
        mode = Models.FixedRho(0.2)
        spec = Models.build_problem_spec(mode)

        params = Models.GapParams(0.5, Models.cached_nodes(8, 4), 0.0)
        cond = spec.conditions(params)
        residual = cond([0.5], [-1.5, -1.5, -2.1, 0.2, 0.2, 1.5, 1.5, 1.5])
        @test length(residual) == 8

        solved = spec.forward_solve(model, 0.5; seed_guess=[-1.5, -1.5, -2.1, 0.2, 0.2, 1.5, 1.5, 1.5], p_num=8, t_num=4)
        @test solved isa NamedTuple
        @test haskey(solved, :converged)
        @test haskey(solved, :residual_norm)
    end
end
