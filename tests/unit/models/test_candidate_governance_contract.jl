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

    attempts = [
        (label=:a, ok=true),
        (label=:b, ok=true),
    ]

    early_stop = Models.execute_attempt_pool(attempts;
        evaluate_attempt=(attempt, idx) -> ((label=attempt.label, idx=idx), attempt.ok),
        on_error=(attempt, idx) -> ((label=attempt.label, idx=idx, err=true), false),
        stop_on_first_success=true,
    )
    @test length(early_stop) == 1
    @test early_stop[1].label == :a

    evaluate_all = Models.execute_attempt_pool(attempts;
        evaluate_attempt=(attempt, idx) -> ((label=attempt.label, idx=idx), attempt.ok),
        on_error=(attempt, idx) -> ((label=attempt.label, idx=idx, err=true), false),
        stop_on_first_success=true,
        evaluate_all_attempts=true,
    )
    @test length(evaluate_all) == 2
    @test evaluate_all[2].label == :b
end
