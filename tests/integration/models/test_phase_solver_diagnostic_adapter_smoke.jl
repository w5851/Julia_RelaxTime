using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const P = Models.pnjl_module()

@testset "phase solver diagnostic adapter smoke" begin
    model = Models.create_model(:PNJL)
    T_fm = 100.0 / 197.327

    @testset "fixed-mu diagnostic readable" begin
        raw = Models.solve_constraint(
            model,
            Models.FixedMu(),
            T_fm;
            μ_fm=0.0,
            fixedmu_use_problem_spec=true,
            diagnostic_level=:summary,
            seed_guess=copy(P.HADRON_SEED_5),
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )
        diag = Models._pm_extract_solver_diagnostic(raw; seed_source=:seed0)
        @test haskey(diag, :attempt_origin)
        @test haskey(diag, :hard_constraint_ok)
        @test haskey(diag, :failed_constraints)
    end

    @testset "trho diagnostic readable" begin
        raw = Models.solve_constraint(
            model,
            Models.FixedRho(0.2),
            T_fm;
            diagnostic_level=:summary,
            seed_guess=copy(P.HADRON_SEED_8),
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            iterations=120,
        )
        diag = Models._pm_extract_solver_diagnostic(raw; seed_source=:seed0)
        status = Models._pm_infer_phase_status(diag)
        @test status in (:valid, :invalid, :unknown)
    end
end
