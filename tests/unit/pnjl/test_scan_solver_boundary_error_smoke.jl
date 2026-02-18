using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))
using .PNJL

@testset "Scan solver boundary error smoke" begin
    tmp_dir = mktempdir()

    # Invalid backend combination: solver_backend=:models requires thermo_backend=:models.
    output_tmu = joinpath(tmp_dir, "tmu_invalid_backend.csv")
    stats_tmu = run_tmu_scan(
        T_values=[90.0],
        mu_values=[10.0],
        xi_values=[0.0],
        output_path=output_tmu,
        overwrite=true,
        resume=false,
        thermo_backend=:legacy,
        solver_backend=:models,
        p_num=8,
        t_num=4,
        iterations=30,
    )

    @test stats_tmu.total == 1
    @test stats_tmu.failure == 1
    lines_tmu = readlines(output_tmu)
    @test length(lines_tmu) == 2
    @test occursin("solver_backend=:models requires thermo_backend=:models", lines_tmu[2])

    output_trho = joinpath(tmp_dir, "trho_invalid_backend.csv")
    stats_trho = run_trho_scan(
        T_values=[90.0],
        rho_values=[0.2],
        xi_values=[0.0],
        output_path=output_trho,
        overwrite=true,
        resume=false,
        thermo_backend=:legacy,
        solver_backend=:models,
        p_num=8,
        t_num=4,
        iterations=30,
    )

    @test stats_trho.total == 1
    @test stats_trho.failure == 1
    lines_trho = readlines(output_trho)
    @test length(lines_trho) == 2
    @test occursin("solver_backend=:models requires thermo_backend=:models", lines_trho[2])
end
