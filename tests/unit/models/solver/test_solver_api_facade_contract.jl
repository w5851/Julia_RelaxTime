using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "solver api facade contract" begin
    T_fm = 100.0 / 197.327
    μ_fm = 0.0

    fixedmu = Models.FixedMu()
    theta_vec = [T_fm, μ_fm]
    theta_named = (T_fm=T_fm, μ_fm=μ_fm)

    @testset "all solver api can be called from Models facade without explicit model" begin
        @test Models.solve(fixedmu, T_fm, μ_fm; p_num=8, t_num=4, residual_norm_max=1e-6) isa Models.SolverResult
        @test Models.solve_multi(fixedmu, T_fm, μ_fm; p_num=8, t_num=4, residual_norm_max=1e-6) isa Models.SolverResult

        raw = Models.solve_constraint(fixedmu, T_fm;
            μ_fm=μ_fm,
            seed_guess=copy(Models.pnjl_module().HADRON_SEED_5),
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )
        @test hasproperty(raw, :converged)

        @test Models.solve_vec(fixedmu, theta_vec; p_num=8, t_num=4, residual_norm_max=1e-6) isa Models.SolverResult
        @test Models.solve_named(fixedmu, theta_named; p_num=8, t_num=4, residual_norm_max=1e-6) isa Models.SolverResult
    end
end
