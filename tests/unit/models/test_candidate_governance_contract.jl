using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "candidate governance normalization helpers" begin
    normalized_ok = Models.normalize_governance_candidate((;
        converged=true,
        pressure=5.0,
        residual_norm=1e-8,
        hard_constraint_ok=true,
        failed_constraints=Symbol[],
    ); seed_index=3)

    @test normalized_ok.seed_index == 3
    @test haskey(normalized_ok, :selection_score)
    @test normalized_ok.converged
    @test normalized_ok.hard_constraint_ok
    @test normalized_ok.failed_constraints == Symbol[]

    normalized_fail = Models.normalize_governance_candidate((;
        converged=false,
        pressure=NaN,
        residual_norm=Inf,
        hard_constraint_ok=false,
        failed_constraints=Symbol[:solver_failed],
    ); seed_index=2)

    @test normalized_fail.seed_index == 2
    @test !normalized_fail.hard_constraint_ok
    @test normalized_fail.failed_constraints == Symbol[:solver_failed]
    @test normalized_fail.pressure == -Inf

    @test_throws ArgumentError Models.normalize_governance_candidate((;
        converged=false,
        pressure=NaN,
        residual_norm=Inf,
        hard_constraint_ok=false,
    ); seed_index=2)

    @test_throws ArgumentError Models.normalize_governance_candidate((;
        converged=true,
        pressure=1.0,
        residual_norm=1e-9,
        hard_constraint_ok=true,
        failed_constraints=Symbol[:residual_too_large],
    ); seed_index=1)

    @test_throws ArgumentError Models.normalize_governance_candidate((;
        converged=false,
        pressure=0.0,
        residual_norm=1e-3,
        hard_constraint_ok=false,
        failed_constraints=Symbol[],
    ); seed_index=1)

    @test !Models.evaluate_candidate_success((; converged=true, residual_norm=1e-2, hard_constraint_ok=false); residual_norm_max=1e-6)
    @test Models.evaluate_candidate_success((; converged=false, residual_norm=1e-8, hard_constraint_ok=true); residual_norm_max=1e-6)
    @test !Models.evaluate_candidate_success((; converged=false, residual_norm=1e-2, hard_constraint_ok=false); residual_norm_max=1e-6)
    @test !Models.evaluate_candidate_success((; converged=true, residual_norm=1e-8, hard_constraint_ok=false); residual_norm_max=1e-6)
    @test Models.classify_attempt_error(ArgumentError("bad")) == :constraint_error
    @test Models.classify_attempt_error(DomainError(1.0, "bad")) == :constraint_error
    @test Models.classify_attempt_error(ErrorException("boom")) == :solver_error
    @test length(Models.normalize_error_message(ErrorException("line1\nline2"))) > 0

    pool = Models.build_seed_pool(Models.FixedEntropy(0.5);
        primary_seed=[1.0, 2.0, 3.0],
        extra_seed_pool=([4.0, 5.0, 6.0], [1.0, 2.0, 3.0]),
        provided_seed_pool=([7.0, 8.0, 9.0],),
        default_seed_pool=([7.0, 8.0, 9.0], [10.0, 11.0, 12.0]),
    )
    @test length(pool) == 4
    @test pool[1].source == :primary
    @test pool[2].source == :extra
    @test pool[3].source == :provided
    @test pool[4].source == :default
end

@testset "candidate governance quality tags and selector contract" begin
    candidates = [
        (
            converged=true,
            pressure=9.0,
            residual_norm=1e-8,
            hard_constraint_ok=true,
            failed_constraints=Symbol[],
            governed_attempt_origin=:primary,
            seed_index=2,
        ),
        (
            converged=true,
            pressure=10.0,
            residual_norm=2e-8,
            hard_constraint_ok=true,
            failed_constraints=Symbol[],
            governed_attempt_origin=:method_rescue,
            seed_index=1,
        ),
        (
            converged=false,
            pressure=11.0,
            residual_norm=1e-3,
            hard_constraint_ok=false,
            failed_constraints=Symbol[:solver_failed],
            seed_index=3,
        ),
    ]

    tagged = Models.normalize_selector_candidates(candidates; residual_norm_max=1e-6)
    @test length(tagged) == 3
    @test tagged[1].quality_tag == :good
    @test tagged[2].quality_tag == :fallback
    @test tagged[3].quality_tag == :bad

    selected = Models.execute_governance_selector(candidates;
        selector=Models.select_pressure_max_candidate,
        residual_norm_max=1e-6,
    )
    @test selected.selected_index == 1
    @test selected.selection_reason == :pressure_max_under_constraints
    @test haskey(selected, :normalized_candidates)
    @test selected.selected_candidate.quality_tag == :good

    @test Models.governance_quality_tag((;
        converged=true,
        residual_norm=1e-8,
        hard_constraint_ok=true,
        governed_attempt_origin=:method_rescue,
    ); residual_norm_max=1e-6) == :fallback

    built = Models.build_governance_candidate((;
        converged=true,
        pressure=8.0,
        residual_norm=1e-8,
        solution=[1.0, 2.0],
        x_state=[1.0, 2.0, 3.0, 0.5, 0.5],
        mu_vec=[0.1, 0.1, 0.1],
        omega=-8.0,
        rho_norm=0.5,
        entropy=0.2,
        energy=1.0,
        masses=[0.3, 0.3, 0.5],
        iterations=1,
    );
        hard_constraint_ok=true,
        failed_constraints=Symbol[],
        seed_index=5,
        residual_norm_max=1e-6,
    )
    @test built.seed_index == 5
    @test built.quality_tag == :good

    selector_checked = Models.execute_governance_selector([
        (
            converged=true,
            pressure=1.0,
            residual_norm=1e-9,
            hard_constraint_ok=true,
            failed_constraints=Symbol[],
        ),
    ]; selector=cands -> begin
        c = cands[1]
        @test hasproperty(c, :seed_index)
        @test hasproperty(c, :selection_score)
        @test hasproperty(c, :quality_tag)
        return (
            selected_index=1,
            selected_candidate=c,
            selection_reason=:pressure_max_under_constraints,
        )
    end)
    @test selector_checked.selected_index == 1
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
        on_error=(attempt, idx, err) -> ((label=attempt.label, idx=idx, err=err !== nothing), false),
        stop_on_first_success=true,
    )
    @test length(early_stop) == 1
    @test early_stop[1].label == :a

    evaluate_all = Models.execute_attempt_pool(attempts;
        evaluate_attempt=(attempt, idx) -> ((label=attempt.label, idx=idx), attempt.ok),
        on_error=(attempt, idx, err) -> ((label=attempt.label, idx=idx, err=err !== nothing), false),
        stop_on_first_success=true,
        evaluate_all_attempts=true,
    )
    @test length(evaluate_all) == 2
    @test evaluate_all[2].label == :b

    errored = Models.execute_attempt_pool([(; label=:x, ok=false)];
        evaluate_attempt=(_, _) -> error("boom"),
        on_error=(attempt, idx, err) -> ((label=attempt.label, idx=idx, err_msg=sprint(showerror, err)), false),
        stop_on_first_success=true,
    )
    @test length(errored) == 1
    @test occursin("boom", errored[1].err_msg)
end
