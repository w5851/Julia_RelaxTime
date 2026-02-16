using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

# Load PNJL module tree (includes scans)
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))
using .PNJL: run_tmu_scan

@testset "TmuScan smoke: single point (solver_backend=:models)" begin
    tmp_dir = mktempdir()

    output = joinpath(tmp_dir, "tmu_scan_solver_backend_models.csv")
    stats = run_tmu_scan(
        T_values=[150.0],
        mu_values=[0.0],
        xi_values=[0.0],
        output_path=output,
        overwrite=true,
        resume=false,
        use_phase_aware=false,
        thermo_backend=:models,
        solver_backend=:models,
        p_num=12,
        t_num=4,
        iterations=80,
    )

    @test isfile(output)
    @test stats.total == 1
    @test stats.success == 1

    lines = readlines(output)
    @test length(lines) == 2
    cols = split(lines[2], ',')
    @test length(cols) >= 4
    @test isfinite(parse(Float64, cols[4]))  # pressure_fm4
end
