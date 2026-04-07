using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const P = Models.pnjl_module()

@testset "solver diagnostic contract" begin
    model = Models.create_model(:PNJL)
    T_fm = 100.0 / 197.327

    @testset "legacy mode parity" begin
        legacy = Models.solve_constraint(
            model,
            Models.FixedMu(),
            T_fm;
            μ_fm=0.0,
            seed_guess=copy(P.HADRON_SEED_5),
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )

        expected_legacy_keys = Set(keys(legacy))
        @test !(:diagnostic in expected_legacy_keys)
        @test Set(keys(legacy)) == expected_legacy_keys
    end

    @testset "diagnostic level none" begin
        mode = Models.FixedEntropy(0.5)
        spec = Models.build_problem_spec(mode)
        result_none = spec.forward_solve(
            model,
            T_fm;
            seed_guess=copy(P.HADRON_SEED_8),
            rho0=0.16,
            diagnostic_level=:none,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            iterations=120,
        )
        @test !haskey(result_none, :diagnostic)
    end

    @testset "diagnostic level summary" begin
        mode = Models.FixedEntropy(0.5)
        spec = Models.build_problem_spec(mode)
        result_summary = spec.forward_solve(
            model,
            T_fm;
            seed_guess=copy(P.HADRON_SEED_8),
            rho0=0.16,
            diagnostic_level=:summary,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            iterations=120,
        )
        @test haskey(result_summary, :diagnostic)
        @test !haskey(result_summary.diagnostic, :candidates)
        @test haskey(result_summary.diagnostic, :attempt_origin)
        @test haskey(result_summary.diagnostic, :seed_source)
        @test haskey(result_summary.diagnostic, :hard_constraint_ok)
        @test haskey(result_summary.diagnostic, :failed_constraints)
        @test haskey(result_summary.diagnostic, :endpoint_cause)
        @test haskey(result_summary.diagnostic, :continuity_distance)
    end

    @testset "diagnostic level full" begin
        mode = Models.FixedRho(0.2)
        spec = Models.build_problem_spec(mode)
        result_full = spec.forward_solve(
            model,
            T_fm;
            seed_guess=copy(P.HADRON_SEED_8),
            diagnostic_level=:full,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            iterations=120,
        )
        @test haskey(result_full, :diagnostic)
        @test haskey(result_full.diagnostic, :candidates)
        @test result_full.diagnostic.candidates isa AbstractVector
    end
end
