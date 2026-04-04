using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const P = Models.pnjl_module()

@testset "firstorder manifold branch stability" begin
    model = Models.create_model(:PNJL)
    mode = Models.FixedEntropy(0.5)
    spec = Models.build_problem_spec(mode)

    T_fm = 150.0 / 197.327
    seed = copy(P.HADRON_SEED_8)

    r1 = spec.forward_solve(
        model,
        T_fm;
        seed_guess=seed,
        rho0=0.16,
        semantic_mode=:constrained_manifold,
        p_num=8,
        t_num=4,
        residual_norm_max=1e-6,
        iterations=120,
    )

    r2 = spec.forward_solve(
        model,
        T_fm;
        seed_guess=seed,
        rho0=0.16,
        semantic_mode=:constrained_manifold,
        p_num=8,
        t_num=4,
        residual_norm_max=1e-6,
        iterations=120,
    )

    @test r1.selection_reason == r2.selection_reason
    @test r1.selected_index == r2.selected_index
    @test isapprox(r1.residual_norm, r2.residual_norm; rtol=1e-10, atol=1e-12)
end
