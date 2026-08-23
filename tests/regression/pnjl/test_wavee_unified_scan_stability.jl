using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Wave-E unified scan stability regression" begin
    @test !isdefined(Models, :scan_workflow_migration_map)
    @test !isdefined(Models, :scan_workflow_migration_status)

    model_kinds = (:PNJL, :NJL, :RPNJL, :Rotation, :GasLiquid)
    for model_kind in model_kinds
        mktempdir() do tmpdir
            tmu_stats = Models.run_tmu_scan(
                T_values=[150.0],
                mu_values=[0.0],
                xi_values=[0.0],
                output_path=joinpath(tmpdir, "tmu.csv"),
                overwrite=true,
                resume=false,
                use_phase_aware=false,
                solver_backend=:models,
                model_kind=model_kind,
                p_num=8,
                t_num=4,
                iterations=20,
            )
            @test tmu_stats.total == 1
            @test tmu_stats.success + tmu_stats.failure == 1

            trho_stats = Models.run_trho_scan(
                T_values=[150.0],
                rho_values=[0.2],
                xi_values=[0.0],
                output_path=joinpath(tmpdir, "trho.csv"),
                overwrite=true,
                resume=false,
                reverse_rho=false,
                seed_policy=:hybrid_continuity,
                constraint_mode=:fixed_rho,
                solver_backend=:models,
                model_kind=model_kind,
                p_num=8,
                t_num=4,
                iterations=20,
            )
            @test trho_stats.total == 1
            @test trho_stats.success + trho_stats.failure == 1
        end
    end

    @testset "magnetic route is explicit" begin
        @test_throws ArgumentError Models.run_tmu_scan(
            T_values=[150.0], mu_values=[0.0], xi_values=[0.0],
            output_path=joinpath(mktempdir(), "magnetic_tmu.csv"),
            model_kind=:PNJLMagnetic,
        )
        @test_throws ArgumentError Models.run_trho_scan(
            T_values=[150.0], rho_values=[0.2], xi_values=[0.0],
            output_path=joinpath(mktempdir(), "magnetic_trho.csv"),
            model_kind=:PNJLMagnetic,
        )
        @test isdefined(Models, :run_magnetic_scan)
    end
end
