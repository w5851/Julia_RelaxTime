using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const P = Models.pnjl_module()

@testset "solver auto backend semantic parity" begin
    mode = P.FixedRho(0.2)
    T_fm = 100.0 / 197.327
    seed = copy(P.HADRON_SEED_8)

    auto = P.TrhoScan._solve_with_models(
        mode,
        T_fm;
        xi=0.0,
        model_kind=:PNJL,
        seed_strategy=P.DefaultSeed(seed, seed, :hadron),
        semantic_mode=:ground_state,
        p_num=8,
        t_num=4,
        residual_norm_max=1e-6,
        iterations=120,
    )

    models = P.TrhoScan._solve_with_models(
        mode,
        T_fm;
        xi=0.0,
        model_kind=:PNJL,
        seed_strategy=P.DefaultSeed(seed, seed, :hadron),
        semantic_mode=:ground_state,
        p_num=8,
        t_num=4,
        residual_norm_max=1e-6,
        iterations=120,
    )

    @test auto.converged == models.converged
    @test isapprox(auto.pressure, models.pressure; rtol=1e-9, atol=1e-10)
    @test isapprox(auto.rho_norm, models.rho_norm; rtol=1e-9, atol=1e-10)
    @test isapprox(auto.residual_norm, models.residual_norm; rtol=1e-9, atol=1e-10)
end
