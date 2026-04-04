using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "auto backend routes to models for PNJL scans" begin
    tmpdir = mktempdir()

    tmu_output = joinpath(tmpdir, "tmu_auto_models.csv")
    trho_output = joinpath(tmpdir, "trho_auto_models.csv")

    tmu_stats = Models.run_tmu_scan(
        T_values=[150.0],
        mu_values=[0.0],
        xi_values=[0.0],
        output_path=tmu_output,
        overwrite=true,
        resume=false,
        solver_backend=:auto,
        auto_pnjl_backend=:models,
        model_kind=:PNJL,
        p_num=8,
        t_num=4,
    )

    trho_stats = Models.run_trho_scan(
        T_values=[150.0],
        rho_values=[0.2],
        xi_values=[0.0],
        output_path=trho_output,
        overwrite=true,
        resume=false,
        reverse_rho=false,
        solver_backend=:auto,
        auto_pnjl_backend=:models,
        model_kind=:PNJL,
        p_num=8,
        t_num=4,
    )

    @test tmu_stats.total == 1
    @test trho_stats.total == 1
    @test isfile(tmu_output)
    @test isfile(trho_output)

    tmu_lines = readlines(tmu_output)
    trho_lines = readlines(trho_output)
    @test length(tmu_lines) == 2
    @test length(trho_lines) == 2
end
