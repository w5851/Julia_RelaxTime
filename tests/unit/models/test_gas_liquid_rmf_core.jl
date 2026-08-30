using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
@testset "GasLiquid RMF core semantics" begin
    p = Models.GasLiquidCoreParams(profile="DiToro_NLrhoDelta")
    @test p.f_sigma ≈ (p.g_sigma / p.m_sigma_inv_fm)^2 atol=1e-12
    @test p.f_omega ≈ (p.g_omega / p.m_omega_inv_fm)^2 atol=1e-12
    @test p.f_rho ≈ (p.g_rho / p.m_rho_inv_fm)^2 atol=1e-12
    @test p.f_delta ≈ (p.g_delta / p.m_delta_inv_fm)^2 atol=1e-12
    @test Models.GasLiquidCoreParams(profile="DiToro_NLrho").f_rho ≈ 0.95
    @test Models.GasLiquidCoreParams(profile="Thesis_NLrho").f_rho ≈ 3.15
    @test Models.GasLiquidCoreParams(profile="Thesis_NLrho").f_delta == 0.0

    T = 15.0
    row_sym = Models.solve_gas_liquid_rmf_point(
        T_MeV=T,
        mode=:fixed_rho,
        rho_B_fm3=0.16,
        alpha=0.0,
        profile="Thesis_NLrho",
        p_num=48,
    )
    @test row_sym.converged
    @test row_sym.rho_B_fm3 ≈ 0.16 atol=2e-8
    @test row_sym.rho_3_fm3 ≈ 0.0 atol=2e-8
    @test row_sym.rho_s3_fm3 ≈ 0.0 atol=2e-8
    @test row_sym.D_inv_fm ≈ 0.0 atol=2e-8
    @test row_sym.M_p_MeV ≈ row_sym.M_n_MeV atol=2e-7
    @test row_sym.entropy_fm3 >= 0.0

    row_delta = Models.solve_gas_liquid_rmf_point(
        T_MeV=T,
        mode=:fixed_rho,
        rho_B_fm3=0.16,
        alpha=0.2,
        profile="DiToro_NLrhoDelta",
        p_num=48,
    )
    @test row_delta.converged
    @test row_delta.rho_B_fm3 ≈ 0.16 atol=2e-8
    @test row_delta.rho_3_fm3 ≈ -0.032 atol=2e-8
    @test abs(row_delta.rho_s3_fm3) > 1e-5
    @test row_delta.D_inv_fm ≈ p.f_delta * row_delta.rho_s3_fm3 atol=2e-7
    @test row_delta.R_inv_fm ≈ p.f_rho * row_delta.rho_3_fm3 atol=2e-8
    # This is the direct consequence of the frozen tau3_p=+1, tau3_n=-1
    # convention and D=f_delta*(rho_s,p-rho_s,n). The opposite mass ordering
    # would require changing that convention and is not silently accepted.
    @test row_delta.M_p_MeV - row_delta.M_n_MeV ≈ -2 * row_delta.D_inv_fm * p.hbarc_MeV_fm atol=2e-5

    thermo = row_sym.thermo
    @test Models.GasLiquidThermodynamics.omega_identity_residual(thermo, T / p.hbarc_MeV_fm) ≈ 0.0 atol=2e-10
end
