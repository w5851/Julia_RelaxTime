using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Production phase pipeline helpers" begin
    @test Models.ProductionPipelineConfig().cep_tol_MeV == 0.1

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
