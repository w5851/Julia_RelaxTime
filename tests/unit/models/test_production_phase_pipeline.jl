using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Production phase pipeline helpers" begin
    @test Models.ProductionPipelineConfig().cep_tol_MeV == 0.1
    @test Models.ProductionPipelineConfig().temperature_resolution_target_MeV == 0.1
    @test Models.RhoHybridVerificationConfig().guard_rule == :extrema_outer_samples_v1
    @test Models.RhoHybridVerificationConfig().targeted_cap == 12
    @test Models.RhoHybridVerificationConfig().candidate_policy == :unique_three_crossing_topology_v1
    @test Models.RhoHybridVerificationConfig().endpoint_policy == :bounded_zero_density_v1

    @testset "nullable geometry scalars remain diagnostic" begin
        unresolved = (
            position_error_MeV=nothing,
            density_error=nothing,
            maxwell_area_gate=nothing,
        )
        normalized = Models._normalize_production_geometry_scalars(unresolved)
        @test normalized.position_error_MeV == Inf
        @test normalized.density_error == Inf
        @test normalized.maxwell_area_gate == Inf
        @test normalized.missing_fields == [
            :position_error_MeV, :density_error, :maxwell_area_gate,
        ]
        @test normalized.normalization_version == :nothing_to_inf_v1

        absent = (
            position_error_MeV=0.01,
            density_error=0.001,
        )
        absent_normalized = Models._normalize_production_geometry_scalars(absent)
        @test absent_normalized.position_error_MeV == 0.01
        @test absent_normalized.density_error == 0.001
        @test absent_normalized.maxwell_area_gate == Inf
        @test absent_normalized.missing_fields == [:maxwell_area_gate]

        finite = (
            position_error_MeV=0.0,
            density_error=0.0025,
            maxwell_area_gate=5e-5,
        )
        finite_normalized = Models._normalize_production_geometry_scalars(finite)
        @test finite_normalized.position_error_MeV == 0.0
        @test finite_normalized.density_error == 0.0025
        @test finite_normalized.maxwell_area_gate == 5e-5
        @test isempty(finite_normalized.missing_fields)

        endpoint = (
            position_error_MeV=0.01,
            density_error=0.001,
            maxwell_area_gate=1e-6,
            rho_hadron=0.0,
            slice_status=:confirmed_first_order,
        )
        endpoint_normalized = Models._normalize_production_geometry_scalars(endpoint)
        @test endpoint_normalized.position_error_MeV == 0.01
        @test endpoint.rho_hadron == 0.0
        @test endpoint.slice_status == :confirmed_first_order

        @test_throws ArgumentError Models._normalize_production_geometry_scalars((
            position_error_MeV="unresolved",
            density_error=0.0,
            maxwell_area_gate=0.0,
        ))
    end

    @testset "hybrid support uses strict outer extrema samples" begin
        sres = Models.SShapeResult(true, 4.0, 2.0, 1.5, 3.5, 2)
        curve = ([0.0, 4.0, 3.0, 2.0, 5.0, 6.0], collect(0.0:1.0:5.0))
        stage_a = (sres=sres,)
        stage_b = (sres=sres, curve=curve, mu_transition=3.0,
            rho_hadron=0.5, rho_quark=4.5)
        support = Models._hybrid_support_grid(collect(0.0:0.05:4.0), stage_a, stage_b)
        @test support.status == :ok
        @test support.low == 0.0
        @test support.high == 4.0
        @test support.mu_low == 2.0
        @test support.mu_high == 4.0
        @test support.guard_rule == :extrema_outer_samples_v1
        @test length(support.grid) <= 12
        @test all(value -> !(value in Set(collect(0.0:0.00625:4.0))), support.grid)
        @test all(0.0 <= value <= 4.0 for value in support.grid)
        @test all(isapprox(value / 0.003125, round(value / 0.003125); atol=1e-8, rtol=0.0) for value in support.grid)
        equal_outer_curve = ([1.0, 1.9, 2.0, 4.0, 3.0, 2.0, 4.1, 4.2],
            [0.0, 0.05, 1.0, 1.5, 2.5, 3.5, 3.55, 4.0])
        equal_outer = Models._hybrid_support_grid(collect(0.0:0.05:4.0), stage_a,
            (sres=sres, curve=equal_outer_curve, mu_transition=3.0))
        @test equal_outer.status == :ok
        @test equal_outer.low == 0.05
        @test equal_outer.high == 3.55
        missing_right = Models._hybrid_support_grid(collect(0.0:0.05:4.0), stage_a,
            (sres=sres, curve=([0.0, 4.0, 3.0, 2.0, 3.5], [0.0, 1.5, 2.5, 3.5, 4.0])))
        @test missing_right.status == :ambiguous_near_critical
        missing_support = Models._hybrid_support_grid(collect(0.0:0.05:4.0),
            (sres=Models.SShapeResult(),),
            (sres=Models.SShapeResult(), curve=nothing))
        @test missing_support.status == :ambiguous_near_critical
    end

    @testset "endpoint-limit route is opt-in and requires a unique public candidate" begin
        sres = Models.SShapeResult(true, 4.0, 2.0, 1.5, 3.5, 2)
        curve = ([2.0, 4.0, 3.0, 2.0, 5.0, 6.0], collect(0.0:1.0:5.0))
        maxwell = Models.MaxwellResult(true, 3.0, 0.001, 4.0, 1e-7, 4,
            Dict{Symbol, Any}(:candidate_count => 1, :crossing_count => 3,
                :endpoint_dependent => true, :crossings => [0.001, 2.0, 4.0]))
        stage_b = (sres=sres, curve=curve, maxwell=maxwell)
        endpoint = Models._hybrid_endpoint_candidate(stage_b, Models.RhoHybridVerificationConfig())
        @test endpoint !== nothing
        @test endpoint.anchor == 0.003125
        no_endpoint = Models._hybrid_endpoint_candidate(stage_b,
            Models.RhoHybridVerificationConfig(endpoint_policy=:other))
        @test no_endpoint === nothing

        v2_curve = ([2.0, 4.0, 3.0, 2.0, 2.5, 4.0], collect(0.0:1.0:5.0))
        v2_sres = Models.SShapeResult(true, 4.0, 2.0, 1.5, 3.5, 2)
        v2_maxwell = Models.MaxwellResult(true, 3.0, 0.001, 4.0, 1e-7, 4,
            Dict{Symbol, Any}(:candidate_count => 1, :crossing_count => 3,
                :endpoint_dependent => true, :crossings => [0.5, 2.0, 4.5]))
        v2_stage_b = (sres=v2_sres, curve=v2_curve, mu_transition=3.0, maxwell=v2_maxwell)
        v2_endpoint = Models._hybrid_endpoint_candidate(v2_stage_b,
            Models.RhoHybridVerificationConfig(endpoint_policy=:three_crossing_endpoint_local_v2))
        @test v2_endpoint !== nothing
        @test v2_endpoint.policy == :three_crossing_endpoint_local_v2
        @test v2_endpoint.anchor == 0.003125
        @test v2_endpoint.left_bracket.low == 0.0
        @test v2_endpoint.left_bracket.high == 1.0
        @test v2_endpoint.right_bracket.low == 4.0
        @test v2_endpoint.right_bracket.high == 5.0

        positive_curve = ([2.0, 4.0, 3.0, 2.0, 2.5, 4.0], collect(0.0:1.0:5.0))
        positive_maxwell = Models.MaxwellResult(true, 3.0, 0.001, 4.0, 1e-7, 4,
            Dict{Symbol, Any}(:candidate_count => 1, :crossing_count => 3,
                :endpoint_dependent => true, :crossings => [0.5, 2.0, 4.5]))
        positive_stage_b = (sres=v2_sres, curve=positive_curve,
            mu_transition=3.0, maxwell=positive_maxwell)
        positive_endpoint = Models._hybrid_endpoint_candidate(positive_stage_b,
            Models.RhoHybridVerificationConfig(endpoint_policy=:three_crossing_endpoint_local_v2))
        @test positive_endpoint !== nothing
        @test positive_endpoint.left_bracket.low == 0.0
        @test positive_endpoint.left_bracket.high == 1.0

        @testset "endpoint support envelope is distinct from active bracket" begin
            envelope = Models._hybrid_endpoint_support_envelope(
                (left_bracket=(low=0.0, high=0.00625),),
                0.003125,
                [0.0015625, 0.00078125, 0.00009765625],
            )
            @test envelope.low == 0.0
            @test envelope.high == 0.00625
            active_bracket = (low=0.00009765625, high=0.0001953125)
            @test active_bracket.low > envelope.low
            @test active_bracket.high < envelope.high
        end
    end

    @testset "endpoint-local route preserves convergence diagnostics" begin
        base_record = (
            axis="rho", level=2, left=1.0, right=2.0, midpoint=2.0,
            position_error_MeV=0.01, density_error=0.001,
            maxwell_area=1e-6, response_rtol=nothing,
            converged=true, reason="ok",
        )
        trace = [
            (
                level=0, rho_midpoint=missing, bracket_low=0.0, bracket_high=0.003125,
                position_error_MeV=0.02, density_error=0.002,
                maxwell_area=2e-6, geometry_converged=true, reason=:anchor,
            ),
            (
                level=1, rho_midpoint=0.0015625, bracket_low=0.0, bracket_high=0.003125,
                position_error_MeV=0.01, density_error=0.001,
                maxwell_area=1e-6, geometry_converged=true, reason="ok",
            ),
        ]
        records = Models._endpoint_local_convergence_records(
            (rho_convergence_records=[base_record],),
            (rho_convergence_records=NamedTuple[],),
            trace,
        )
        @test length(records) == 3
        @test records[end - 1].midpoint == 0.0015625
        @test records[end].midpoint == 0.0015625
        @test records[end].axis == "rho"
        @test records[end].converged
    end

    @testset "hybrid policy validation is explicit and opt-in" begin
        @test_throws ArgumentError Models.run_production_phase_pipeline(
            :PNJL; T_start=150.0, T_end=150.0, dT=1.0, rho_grid=[0.0, 0.05],
            rho_refinement_policy=:rho_support_hybrid, cep_max_refine_level_rho=1,
            rho_geometry_convergence=true, promote_reference=false,
        )
        @test_throws ArgumentError Models.run_production_phase_pipeline(
            :PNJL; T_start=150.0, T_end=150.0, dT=1.0, rho_grid=[0.0, 0.05],
            rho_refinement_policy=:rho_support_hybrid, cep_max_refine_level_rho=4,
            rho_geometry_convergence=false, promote_reference=false,
        )
        @test_throws ArgumentError Models.run_production_phase_pipeline(
            :PNJL; T_start=150.0, T_end=150.0, dT=1.0, rho_grid=[0.0, 0.05],
            rho_refinement_policy=:rho_support_hybrid, cep_max_refine_level_rho=4,
            rho_geometry_convergence=true, rho_support_targeted_cap=13,
            promote_reference=false,
        )
        @test_throws ArgumentError Models.run_production_phase_pipeline(
            :PNJL; T_start=150.0, T_end=150.0, dT=1.0, rho_grid=[0.0, 0.05],
            rho_refinement_policy=:rho_support_hybrid, cep_max_refine_level_rho=4,
            rho_geometry_convergence=true,
            rho_hybrid_verification=Models.RhoHybridVerificationConfig(point_ranking_version=:other),
            promote_reference=false,
        )
    end

    temps = Models._production_temperature_grid(150.0, 160.0, 5.0)
    @test temps == [150.0, 155.0, 160.0]

    temps_tail = Models._production_temperature_grid(150.0, 159.0, 5.0)
    @test temps_tail == [150.0, 155.0, 159.0]

    records = [
        (T_MeV=150.0, status=:valid, mu_transition_MeV=250.0, area_residual=1e-5, reason="ok", level=0),
        (T_MeV=155.0, status=:unknown, mu_transition_MeV=245.0, area_residual=3e-4, reason="gray", level=1),
        (T_MeV=160.0, status=:invalid, mu_transition_MeV=NaN, area_residual=NaN, reason="no_s_shape", level=0),
    ]
    sweep = Models._materialize_sweep_result(records, true, 155.0, 1, 1)
    @test sweep.first_point_fallback
    @test sweep.fallback_start_T_MeV == 155.0
    @test sweep.unknown_count == 1
    @test sweep.forced_invalid_count == 1
    @test sweep.statuses == [:valid, :unknown, :invalid]

    bracket = Models._find_production_bracket(records)
    @test bracket !== nothing
    @test bracket.T_low == 150.0
    @test bracket.T_high == 160.0

    @testset "unknown budget stops frontier refinement without relabelling" begin
        cfg_budget = Models.ProductionPipelineConfig(
            temperature_resolution_target_MeV=0.5,
            cep_max_bisect_iter=8,
            unknown_budget=0,
        )
        calls = Ref(0)
        frontiers = Models._refine_production_cep_frontiers(
            (T_low=130.0, mu_low=295.0, T_high=140.0),
            T -> begin
                calls[] += 1
                (slice_status=:ambiguous, mu_transition=nothing)
            end,
            cfg_budget;
            unknown_count_fn=() -> 1,
        )
        @test frontiers.budget_exhausted
        @test frontiers.low.T == 130.0
        @test frontiers.high.T == 140.0
        @test calls[] == 0
    end

    cfg = Models.ProductionPipelineConfig(T_start=120.0, T_end=150.0, dT_initial=5.0)
    cep = Models.CEPResult(found=true, T_cep_MeV=126.09375, uncertainty_T_MeV=0.3125)
    boundary = [
        (T_MeV=120.0, mu_transition_MeV=290.0, rho_hadron=1.0, rho_quark=2.0, area_residual=1e-5, converged=true),
        (T_MeV=125.0, mu_transition_MeV=299.7, rho_hadron=1.1, rho_quark=2.1, area_residual=1e-5, converged=true),
    ]
    crossover_bounds = Models._production_crossover_temperature_bounds(cfg, cep, boundary)
    @test crossover_bounds.T_min_MeV == 125.0
    @test crossover_bounds.T_max_MeV == 150.0

    capped_cfg = Models.ProductionPipelineConfig(
        T_start=120.0,
        T_end=260.0,
        dT_initial=5.0,
        crossover_T_max_MeV=240.0,
    )
    capped_bounds = Models._production_crossover_temperature_bounds(capped_cfg, cep, boundary)
    @test capped_bounds.T_max_MeV == 240.0

    @testset "temperature interpolation refinement preserves sorted evaluation" begin
        calls = Float64[]
        evaluator = function (T::Float64)
            push!(calls, T)
            shift = 0.2 * (T - 1.0)^2
            sres = Models.SShapeResult(true, 310.0 + shift, 290.0 + shift, 0.8, 2.2, 2)
            return (
                status=:valid,
                reason="ok",
                mu_transition=300.0 + shift,
                rho_hadron=1.0,
                rho_quark=2.0,
                area_residual=5e-5,
                sres=sres,
            )
        end
        adaptive_cfg = Models.ProductionPipelineConfig(
            adaptive_temperature=true,
            temperature_max_refine_level=1,
            temperature_position_tol_MeV=0.05,
            temperature_density_tol=0.01,
            temperature_maxwell_area_tol=1e-4,
        )
        resolved, audit = Models._adaptive_production_temperature_refinement!(
            [0.0, 2.0], evaluator, adaptive_cfg,
        )
        @test resolved == [0.0, 1.0, 2.0]
        @test length(audit) == 1
        @test !audit[1].converged
        @test audit[1].position_error_MeV > adaptive_cfg.temperature_position_tol_MeV
        @test calls[1:2] == [0.0, 2.0]
        @test calls[3] == 1.0
    end
end
