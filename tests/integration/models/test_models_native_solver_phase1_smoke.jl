using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const HBARC_MEV_FM = 197.327

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Models native solver phase1 smoke" begin
    @test Models.HADRON_SEED_5 isa AbstractVector
    @test Models.MultiSeed() isa Models.SeedStrategy
    @test Models.PhaseAwareSeed() isa Models.SeedStrategy

    model = Models.create_model(:PNJL)
    T_fm = 100.0 / HBARC_MEV_FM
    μ_fm = 20.0 / HBARC_MEV_FM

    result = Models.solve(model, Models.FixedMu(), T_fm, μ_fm;
        seed_strategy=Models.DefaultSeed(),
        p_num=8,
        t_num=4,
    )

    @test result isa Models.SolverResult
    @test result.converged
    @test isfinite(result.pressure)

    multi_result = Models.solve_multi(model, Models.FixedMu(), T_fm, μ_fm;
        seed_strategy=Models.MultiSeed(),
        p_num=8,
        t_num=4,
    )

    @test multi_result isa Models.SolverResult
    @test multi_result.converged

    deriv_result = Models.solve_with_derivatives(T_fm, μ_fm;
        order=1,
        p_num=8,
        t_num=4,
    )
    @test haskey(deriv_result, :x)
    @test haskey(deriv_result, :dx_dT)
    @test haskey(deriv_result, :dx_dμ)
end
