using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "candidate governance contract" begin
    @test isdefined(Models, :build_candidate_context)
    @test isdefined(Models, :evaluate_hard_constraints)
    @test isdefined(Models, :select_pressure_max_candidate)

    context = Models.build_candidate_context(
        Models.FixedMu();
        continuity_seed_available=true,
        phase_hint=:hadron,
        residual_norm_max=1e-6,
        prefer_continuity=true,
    )

    @test context.mode isa Models.FixedMu
    @test context.continuity_seed_available
    @test context.phase_hint == :hadron
    @test context.residual_norm_max ≈ 1e-6
    @test context.prefer_continuity

    rules = [
        (c, p, ctx) -> (c.residual_norm <= p.max_residual, :residual_too_large),
        (c, p, ctx) -> (!ctx.prefer_continuity || c.seed_source == :continuity, :continuity_preferred),
    ]

    params = (max_residual=1e-4,)
    candidate_ok = (residual_norm=1e-6, seed_source=:continuity)
    candidate_bad = (residual_norm=1e-3, seed_source=:default)

    pass_ok, failed_ok = Models.evaluate_hard_constraints(candidate_ok, rules, params, context)
    @test pass_ok
    @test isempty(failed_ok)

    pass_bad, failed_bad = Models.evaluate_hard_constraints(candidate_bad, rules, params, context)
    @test !pass_bad
    @test :residual_too_large in failed_bad
    @test :continuity_preferred in failed_bad

    candidates = [
        (pressure=9.0, residual_norm=1e-7, hard_constraint_ok=true, failed_constraints=Symbol[], converged=true),
        (pressure=10.0, residual_norm=1e-8, hard_constraint_ok=false, failed_constraints=[:residual_too_large], converged=true),
    ]

    selected = Models.select_pressure_max_candidate(candidates, params, context)
    @test selected.selected_index == 1
    @test selected.selection_reason == :pressure_max_under_constraints
end
