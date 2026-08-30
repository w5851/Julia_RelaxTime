using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "GasLiquid RMF scans and contracts" begin
    point = Models.solve_gas_liquid_rmf_point(T_MeV=15.0, mu_B_MeV=700.0, p_num=24, run_id="test-point", point_id="p1")
    rebuilt = Models.build_gas_liquid_result_row(point)
    @test rebuilt.schema_version == "gas_liquid_rmf_row_v1"
    @test rebuilt.point_id == "p1"
    @test rebuilt.mu_p_MeV ≈ point.mu_p_MeV atol=1e-10

    tmu = Models.run_gas_liquid_tmu_scan([15.0, 20.0], [700.0]; p_num=24, run_id="test-tmu")
    @test tmu.schema_version == "gas_liquid_rmf_row_v1"
    @test length(tmu.rows) == 2
    @test tmu.manifest.row_count == 2
    @test all(haskey(row, :rho_p_fm3) && haskey(row, :rho_n_fm3) for row in tmu.rows)
    @test all(!haskey(row, :state) && !haskey(row, :model) && !haskey(row, :result) for row in tmu.rows)
    @test all(row.formal_status == :diagnostic_only for row in tmu.rows)
    @test tmu.manifest.formal_status == :diagnostic_only

    trho = Models.run_gas_liquid_trho_scan([15.0], [0.08, 0.16]; alpha=0.2, p_num=24, run_id="test-trho")
    @test length(trho.rows) == 2
    @test all(row.mode == :fixed_rho for row in trho.rows)
    @test all((isapprox(row.rho_B_fm3, row.rho_p_fm3 + row.rho_n_fm3; atol=1e-7) for row in trho.rows if row.converged))
    @test all(!haskey(row, :state) && !haskey(row, :model) && !haskey(row, :result) for row in trho.rows)
    @test all(row.schema_version == "gas_liquid_rmf_row_v1" for row in trho.rows)
    @test haskey(trho.manifest, :source_matrix_ids)
    @test haskey(trho.manifest, :failed_points)
end
