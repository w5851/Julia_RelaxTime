using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

# Load Models unified scan entrypoint
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
using .Models: run_trho_scan

@testset "TrhoScan smoke: single point (legacy/models)" begin
    tmp_dir = mktempdir()

    for backend in (:legacy, :models)
        output = joinpath(tmp_dir, "trho_scan_$(backend).csv")
        stats = run_trho_scan(
            T_values=[150.0],
            rho_values=[0.2],
            xi_values=[0.0],
            output_path=output,
            overwrite=true,
            resume=false,
            reverse_rho=false,
            thermo_backend=backend,
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
end
