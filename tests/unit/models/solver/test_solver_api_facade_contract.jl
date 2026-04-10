using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "solver api facade contract" begin
    T_fm = 100.0 / 197.327
    μ_fm = 0.0
    model = Models.create_model(:PNJL)

    fixedmu = Models.FixedMu()
    theta_vec = [T_fm, μ_fm]
    theta_named = (T_fm=T_fm, μ_fm=μ_fm)

    @testset "solver api facade entrypoints stay callable" begin
        @test Models.solve(fixedmu, T_fm, μ_fm; p_num=8, t_num=4, residual_norm_max=1e-6) isa Models.SolverResult
        @test Models.solve_multi(fixedmu, T_fm, μ_fm; p_num=8, t_num=4, residual_norm_max=1e-6) isa Models.SolverResult

        raw = Models.solve_constraint(model, fixedmu, T_fm;
            μ_fm=μ_fm,
            seed_guess=copy(Models.pnjl_module().HADRON_SEED_5),
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )
        @test hasproperty(raw, :converged)
        @test hasproperty(raw, :solution)
        @test hasproperty(raw, :residual_norm)

        @test Models.solve_vec(model, fixedmu, theta_vec; p_num=8, t_num=4, residual_norm_max=1e-6) isa Models.SolverResult
        @test Models.solve_named(model, fixedmu, theta_named; p_num=8, t_num=4, residual_norm_max=1e-6) isa Models.SolverResult
    end
end
