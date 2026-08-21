using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
@testset "GasLiquid symmetric RMF fixed point" begin
    row = Models.solve_gas_liquid_rmf_point(
        T_MeV=15.0,
        mode=:fixed_rho,
        rho_B_fm3=0.16,
        alpha=0.0,
        profile="Thesis_NLrho",
        p_num=64,
        p_max_inv_fm=10.0,
    )
    @test row.converged
    @test row.solver_status == :converged
    @test row.rho_B_fm3 ≈ 0.1600000000031738 atol=2e-8
    @test row.M_p_MeV ≈ 706.6501287579312 rtol=3e-7
    @test row.M_n_MeV ≈ 706.6501287579312 rtol=3e-7
    @test row.S_inv_fm ≈ 1.1774865797422849 rtol=3e-7
    @test row.W_inv_fm ≈ 0.8676800000172115 rtol=3e-7
    @test row.pressure_fm4 ≈ 0.14117735331239184 rtol=3e-7
    @test row.entropy_fm3 ≈ 0.2300055068202937 rtol=3e-7
    @test row.energy_fm4 ≈ 0.6223691942336058 rtol=3e-7
    @test row.formal_status == :diagnostic_only
end
