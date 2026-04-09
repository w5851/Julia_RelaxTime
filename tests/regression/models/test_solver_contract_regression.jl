using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const HBARC_MEV_FM = 197.327

@testset "solver contract regression" begin
    model = Models.create_model(:PNJL)

    fixedmu_raw = Models.solve_constraint(
        model,
        Models.FixedMu(),
        120.0 / HBARC_MEV_FM;
        μ_fm=0.0,
        seed_guess=copy(Models.pnjl_module().HADRON_SEED_5),
        p_num=8,
        t_num=4,
        residual_norm_max=1e-6,
    )
    fixedmu_res = Models.coerce_solver_result(Models.FixedMu(), fixedmu_raw)
    fixedmu_view = Models.solver_result_view(fixedmu_res)

    @test fixedmu_view.contract_version == Models.SOLVER_CONTRACT_VERSION_V1
    for key in Models.SOLVER_RESULT_REQUIRED_FIELDS
        @test haskey(fixedmu_view, key)
    end

    fixedrho_mode = Models.FixedRho(0.2)
    spec = Models.build_problem_spec(fixedrho_mode)
    fixedrho_raw = spec.forward_solve(
        model,
        120.0 / HBARC_MEV_FM;
        seed_guess=copy(Models.pnjl_module().HADRON_SEED_8),
        p_num=8,
        t_num=4,
        residual_norm_max=1e-6,
        diagnostic_level=:summary,
    )

    summary = Models.coerce_solver_diagnostic_summary(fixedrho_raw.diagnostic)
    public_diag = Models.to_public_namedtuple(summary)

    @test Models.solver_diagnostic_version(summary) == Models.SOLVER_DIAGNOSTIC_VERSION_V1
    @test public_diag.diagnostic_version == Models.SOLVER_DIAGNOSTIC_VERSION_V1
    @test !haskey(public_diag, :selection_reason_source)
    @test haskey(fixedrho_raw.diagnostic, :selection_reason_source)
end
