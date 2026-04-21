using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const SUPPORTED_MODEL_KINDS = (:PNJL, :NJL, :RPNJL, :PNJLMagnetic, :Rotation, :GasLiquid)

@testset "Wave-E model_kind x scan_mode matrix smoke" begin
    @testset "supported model kinds run in both scan modes" begin
        for model_kind in SUPPORTED_MODEL_KINDS
            tmpdir = mktempdir()

            tmu_output = joinpath(tmpdir, "tmu_$(model_kind).csv")
            tmu_stats = Models.run_tmu_scan(
                T_values=[150.0],
                mu_values=[0.0],
                xi_values=[0.0],
                output_path=tmu_output,
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
            @test isfile(tmu_output)

            trho_output = joinpath(tmpdir, "trho_$(model_kind).csv")
            trho_stats = Models.run_trho_scan(
                T_values=[150.0],
                rho_values=[0.2],
                xi_values=[0.0],
                output_path=trho_output,
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
            @test isfile(trho_output)
        end
    end

    @testset "pnjl_aniso stays parameterized mode" begin
        err_tmu = try
            Models.run_tmu_scan(
                T_values=[150.0],
                mu_values=[0.0],
                xi_values=[0.1],
                output_path=joinpath(mktempdir(), "aniso.csv"),
                overwrite=true,
                resume=false,
                use_phase_aware=false,
                solver_backend=:models,
                model_kind=:pnjl_aniso,
                p_num=8,
                t_num=4,
                iterations=20,
            )
            nothing
        catch exc
            exc
        end

        @test err_tmu isa ArgumentError
        msg_tmu = lowercase(sprint(showerror, err_tmu))
        @test occursin("pnjl_aniso", msg_tmu)
        @test occursin("profile", msg_tmu)

        err_trho = try
            Models.run_trho_scan(
                T_values=[150.0],
                rho_values=[0.2],
                xi_values=[0.1],
                output_path=joinpath(mktempdir(), "aniso_trho.csv"),
                overwrite=true,
                resume=false,
                reverse_rho=false,
                seed_policy=:hybrid_continuity,
                constraint_mode=:fixed_rho,
                solver_backend=:models,
                model_kind=:pnjl_aniso,
                p_num=8,
                t_num=4,
                iterations=20,
            )
            nothing
        catch exc
            exc
        end

        @test err_trho isa ArgumentError
        msg_trho = lowercase(sprint(showerror, err_trho))
        @test occursin("pnjl_aniso", msg_trho)
        @test occursin("profile", msg_trho)
    end
end
