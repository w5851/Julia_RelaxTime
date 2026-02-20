using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

# Load PNJL module tree (includes scans)
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))
using .PNJL: run_trho_scan

@testset "TrhoScan smoke: single point (solver_backend=:models)" begin
    tmp_dir = mktempdir()

    output = joinpath(tmp_dir, "trho_scan_solver_backend_models.csv")
    stats = run_trho_scan(
        T_values=[150.0],
        rho_values=[0.2],
        xi_values=[0.0],
        output_path=output,
        overwrite=true,
        resume=false,
        reverse_rho=false,
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
    header = split(lines[1], ',')
    cols = split(lines[2], ',')
    @test length(cols) >= 8

    idx_mu_u = findfirst(==("mu_u_MeV"), header)
    idx_pressure = findfirst(==("pressure_fm4"), header)
    @test idx_mu_u !== nothing
    @test idx_pressure !== nothing

    @test isfinite(parse(Float64, cols[idx_mu_u]))
    @test isfinite(parse(Float64, cols[idx_pressure]))
end

@testset "TrhoScan smoke: asymmetric mode (solver_backend=:models)" begin
    tmp_dir = mktempdir()

    for kind in (:PNJL, :RPNJL)
        output = joinpath(tmp_dir, "trho_scan_solver_backend_models_asym_$(kind).csv")
        stats = run_trho_scan(
            T_values=[150.0],
            rho_values=[0.05],
            xi_values=[0.0],
            output_path=output,
            overwrite=true,
            resume=false,
            reverse_rho=false,
            constraint_mode=:fixed_asymmetric_rho,
            asym_ud_ratio_target=0.876,
            asym_s_target=0.0,
            thermo_backend=:models,
            solver_backend=:models,
            model_kind=kind,
            p_num=12,
            t_num=4,
            iterations=120,
        )

        @test isfile(output)
        @test stats.total == 1
        @test stats.success == 1
    end
end

@testset "TrhoScan smoke: single point (solver_backend=:auto, thermo_backend=:models)" begin
    tmp_dir = mktempdir()

    output = joinpath(tmp_dir, "trho_scan_solver_backend_auto.csv")
    stats = run_trho_scan(
        T_values=[150.0],
        rho_values=[0.2],
        xi_values=[0.0],
        output_path=output,
        overwrite=true,
        resume=false,
        reverse_rho=false,
        thermo_backend=:models,
        solver_backend=:auto,
        p_num=12,
        t_num=4,
        iterations=80,
    )

    @test isfile(output)
    @test stats.total == 1
    @test stats.success == 1
end
