using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "problem spec fixedrho forward_solve smoke" begin
    model = Models.create_model(:PNJL)
    mode = Models.FixedRho(0.2)
    spec = Models.build_problem_spec(mode)

    result = Models.solve_constraint(
        model,
        mode,
        90.0 / 197.327;
        problem_spec=spec,
        seed_guess=Models.pnjl_module().HADRON_SEED_8,
        p_num=8,
        t_num=4,
        residual_norm_max=1e-6,
    )

    @test result isa NamedTuple
    @test haskey(result, :converged)
    @test haskey(result, :rho_norm)
    @test haskey(result, :residual_norm)
end
