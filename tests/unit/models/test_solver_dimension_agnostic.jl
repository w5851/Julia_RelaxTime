using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const P = Models.pnjl_module()

@testset "solver dimension-agnostic contract" begin
    @testset "legacy compat APIs removed" begin
        @test !isdefined(Models, :solve_fixedmu_constraint)
        @test !isdefined(Models, :solve_fixedrho_constraint)
        @test !isdefined(Models, :solve_fixedentropy_constraint)
        @test !isdefined(Models, :solve_fixedsigma_constraint)
        @test !isdefined(Models, :solve_fixedasymrho_constraint)
        @test !isdefined(Models, :solver_migration_map)
        @test !isdefined(Models, :solver_migration_status)
        @test !isdefined(Models, :legacy_transport_c_p)
        @test isdefined(Models, :solve_constraint)
        @test :create_implicit_solver ∉ names(Models)
        @test :solve_with_root_diagnostics ∉ names(Models)
        @test :solve_weighted_block_fallback ∉ names(Models)
    end

    @testset "SolverResult accepts non-5/3 state and mu" begin
        result = Models.SolverResult(
            Models.FixedMu(),
            true,
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            [0.1, 0.2, 0.3, 0.4],
            [0.01, 0.02],
            -1.0,
            1.0,
            0.0,
            0.0,
            0.0,
            [0.31, 0.32, 0.33],
            3,
            1e-8,
            0.0,
        )

        @test length(result.x_state) == 4
        @test length(result.mu_vec) == 2
    end

    @testset "build_conditions can dispatch with schema" begin
        params = P.GapParams(0.5, P.cached_nodes(24, 6), 0.0)
        mode = P.FixedRho(1.0)
        schema = Models.schema_for_model(:PNJL)

        θ = [0.5]
        x = [-1.5, -1.5, -2.1, 0.2, 0.2, 1.5, 1.5, 1.5]

        legacy = P.build_conditions(mode, params)
        schema_driven = P.build_conditions(mode, params, schema; mu_dim=3)

        @test schema_driven(θ, x) ≈ legacy(θ, x) rtol=1e-12 atol=1e-12
    end

    @testset "solver backend defaults are model-agnostic" begin
        mktempdir() do tmpdir
            out_tmu = joinpath(tmpdir, "_tmp_models_default_tmu.csv")
            out_trho = joinpath(tmpdir, "_tmp_models_default_trho.csv")

            stats_tmu = Models.run_tmu_scan(
                T_values=[150.0],
                mu_values=[0.0],
                xi_values=[0.0],
                model_kind=:RPNJL,
                output_path=out_tmu,
                overwrite=true,
                resume=false,
                p_num=8,
                t_num=4,
            )
            @test stats_tmu.total == 1
            @test stats_tmu.success + stats_tmu.failure == 1

            stats_trho = Models.run_trho_scan(
                T_values=[150.0],
                rho_values=[0.2],
                xi_values=[0.0],
                model_kind=:RPNJL,
                output_path=out_trho,
                overwrite=true,
                resume=false,
                p_num=8,
                t_num=4,
            )
            @test stats_trho.total == 1
            @test stats_trho.success + stats_trho.failure == 1
        end
    end
end
