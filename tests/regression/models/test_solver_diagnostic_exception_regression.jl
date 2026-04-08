using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "solver diagnostic exception regression" begin
    model = Models.create_model(:PNJL)
    mode = Models.FixedEntropy(0.5)
    spec = Models.build_problem_spec(mode)
    T_fm = 100.0 / 197.327
    seed = copy(Models.pnjl_module().HADRON_SEED_8)

    result = spec.forward_solve(
        model,
        T_fm;
        seed_guess=seed,
        p_num=8,
        t_num=4,
        residual_norm_max=1e-6,
        iterations=80,
        diagnostic_level=:summary,
    )

    @test haskey(result, :diagnostic)
    @test result.diagnostic.error_kind == :constraint_error
    @test occursin("rho0 is required", result.diagnostic.error_msg)
    @test result.selection_reason in (:pressure_max_under_constraints, :residual_min_under_constraints, :no_candidate_passed_constraints)
end
