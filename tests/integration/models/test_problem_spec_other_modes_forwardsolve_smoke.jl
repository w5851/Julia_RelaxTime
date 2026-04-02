using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const HBARC_MEV_FM = 197.327

@testset "problem spec forward_solve smoke for entropy/sigma/asymrho" begin
    model = Models.create_model(:PNJL)
    T_fm = 100.0 / HBARC_MEV_FM

    modes = (
        Models.FixedEntropy(0.5),
        Models.FixedSigma(10.0),
        Models.FixedAsymmetricRho(0.05, 1.0, 0.0),
    )

    for mode in modes
        spec = Models.ProblemSpec(
            mode;
            x_dim=Models.state_dim(mode),
            theta_dim=Models.param_dim(mode),
            forward_solve=(m, t; fwd_kwargs...) -> (
                converged=true,
                residual_norm=0.0,
                used_problem_spec=true,
                mode_tag=string(typeof(mode)),
            ),
        )
        result = Models.solve_constraint(
            model,
            mode,
            T_fm;
            problem_spec=spec,
        )

        @test result isa NamedTuple
        @test result.converged
        @test result.used_problem_spec
        @test result.residual_norm == 0.0
        @test occursin("Fixed", result.mode_tag)
    end
end
