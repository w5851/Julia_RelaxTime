using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
const GUARD_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "analysis", "pnjl_cep_hybrid_extrema_guard_feasibility.jl")
if !isdefined(Main, :StageCExtremaGuardFeasibility)
    include(GUARD_SCRIPT)
end
const GUARD = Main.StageCExtremaGuardFeasibility

function _synthetic_curve(mu)
    points = [GUARD.Point(0.0, 60.0, Float64(i - 1), Float64(value), 0.0)
        for (i, value) in enumerate(mu)]
    (rho=Float64.(0:(length(mu) - 1)), mu=Float64.(mu), points=points)
end

function _synthetic_eval(curve; mu_max=4.0, mu_min=2.0, rho_max=1.5, rho_min=3.5)
    (
        status=:valid, reason="ok", mu_transition=3.0,
        rho_hadron=0.5, rho_quark=4.5, area_residual=1e-6,
        sres=Models.SShapeResult(true, mu_max, mu_min, rho_max, rho_min, 2),
        maxwell=Models.MaxwellResult(), curve=curve,
    )
end

@testset "Stage-C extrema outer-sample guard feasibility" begin
    @testset "strict first outer samples" begin
        curve = _synthetic_curve([0.0, 4.0, 3.0, 2.0, 5.0, 6.0])
        guard = GUARD._extrema_outer_guard(curve, _synthetic_eval(curve))
        @test guard.status == :ok
        @test guard.guard_low == 0.0
        @test guard.guard_high == 4.0
        @test guard.guard_low_mu < guard.mu_low
        @test guard.guard_high_mu > guard.mu_high
    end

    @testset "equal extrema samples are skipped" begin
        curve = _synthetic_curve([2.0, 4.0, 3.0, 2.0, 4.0, 5.0])
        guard = GUARD._extrema_outer_guard(curve, _synthetic_eval(curve))
        @test guard.status == :ambiguous_near_critical
        @test guard.reason == "missing_left_strict_outer_sample"
    end

    @testset "missing outer side and multiple topology remain ambiguous" begin
        curve = _synthetic_curve([0.0, 4.0, 3.0, 2.0])
        guard = GUARD._extrema_outer_guard(curve, _synthetic_eval(curve))
        @test guard.status == :ambiguous_near_critical

        multi = _synthetic_curve([0.0, 4.0, 3.0, 2.0, 5.0, 4.0, 6.0, 5.0, 7.0])
        multi_guard = GUARD._extrema_outer_guard(multi, _synthetic_eval(multi))
        @test multi_guard.reason == "unstable_or_multiple_s_topology"
    end

    @testset "Stage-B curve is retained in the Stage-C union" begin
        base = _synthetic_curve([0.0, 4.0, 3.0, 2.0, 5.0, 6.0])
        extra = [GUARD.Point(0.0, 60.0, 0.5 + i / 100, 2.0 + i / 100, 0.0) for i in 1:20]
        merged = GUARD._curve(vcat(base.points, extra))
        @test merged !== nothing
        @test length(merged.points) == length(base.points) + length(extra)
        @test all(point -> point in merged.points, base.points)
    end

    @testset "target selection is deterministic and capped" begin
        curve = _synthetic_curve([0.0, 4.0, 3.0, 2.0, 5.0, 6.0])
        evaluation = _synthetic_eval(curve)
        guard = GUARD._extrema_outer_guard(curve, evaluation)
        pool = [GUARD.Point(0.0, 60.0, i / 100, 2.0 + i / 1000, 0.0) for i in 1:400]
        first = GUARD._select_stage_c_points(pool, curve, evaluation, guard, 12)
        second = GUARD._select_stage_c_points(pool, curve, evaluation, guard, 12)
        @test first == second
        @test length(first) == 12
        @test length(unique(point.rho for point in first)) == 12
    end

    @testset "default production policies are unchanged" begin
        @test Models.ProductionPipelineConfig().rho_refinement_policy == :uniform_nested
        @test Models.ProductionPipelineConfig().rho_support_targeted_cap == 12
        @test GUARD.STAGE_B_FINE == 0.00625
        @test GUARD.STAGE_C_FINE == 0.003125
        @test GUARD.COMPARISON_EPS == 32eps(Float64)
    end
end
