using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const HBARC_MEV_FM = 197.327

@testset "problem spec forward_solve smoke for entropy/sigma/asymrho" begin
    model = Models.create_model(:PNJL)
    T_fm = 100.0 / HBARC_MEV_FM
    seed = copy(Models.pnjl_module().HADRON_SEED_8)

    modes = (
        Models.FixedEntropy(0.5),
        Models.FixedSigma(10.0),
        Models.FixedAsymmetricRho(0.05, 1.0, 0.0),
    )

    for mode in modes
        spec = Models.build_problem_spec(mode)
        result = Models.solve_constraint(
            model,
            mode,
            T_fm;
            problem_spec=spec,
            seed_guess=seed,
            rho0=0.16,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            iterations=120,
        )

        @test result isa NamedTuple
        @test haskey(result, :converged)
        @test haskey(result, :residual_norm)
        @test haskey(result, :selection_reason)
        @test haskey(result, :candidate_count)
        @test isfinite(result.residual_norm)
        @test result.selection_reason in (:pressure_max_under_constraints, :no_candidate_passed_constraints)
        @test result.candidate_count >= 1
    end
end
