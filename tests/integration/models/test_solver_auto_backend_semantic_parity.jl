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

@testset "solver fixedmu problemspec parity guard" begin
    model = Models.create_model(:PNJL)
    mode = Models.FixedMu()
    T_fm = 100.0 / 197.327
    μ_fm = 0.0
    seed = copy(P.HADRON_SEED_5)

    legacy = Models.solve_constraint(
        model,
        mode,
        T_fm;
        μ_fm=μ_fm,
        seed_guess=seed,
        p_num=8,
        t_num=4,
        residual_norm_max=1e-6,
        fixedmu_use_problem_spec=false,
    )

    spec_path = Models.solve_constraint(
        model,
        mode,
        T_fm;
        μ_fm=μ_fm,
        seed_guess=seed,
        p_num=8,
        t_num=4,
        residual_norm_max=1e-6,
        fixedmu_use_problem_spec=true,
    )

    @test haskey(legacy, :fixedmu_problem_spec_active)
    @test haskey(spec_path, :fixedmu_problem_spec_active)
    @test legacy.fixedmu_problem_spec_active === false
    @test spec_path.fixedmu_problem_spec_active === true

    @test legacy.converged == spec_path.converged
    @test isapprox(legacy.pressure, spec_path.pressure; rtol=1e-9, atol=1e-10)
    @test isapprox(legacy.rho_norm, spec_path.rho_norm; rtol=1e-9, atol=1e-10)
    @test isapprox(legacy.residual_norm, spec_path.residual_norm; rtol=1e-9, atol=1e-10)
end
