using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "GasLiquid RMF formula mapping" begin
    p = Models.GasLiquidCoreParams(profile="DiToro_NLrhoDelta")
    st = Models.GasLiquidEquationSet.GasLiquidState(0.2, -0.04, 3.2, 3.0)
    dens = Models.GasLiquidEquationSet.density_bundle(st, 15.0 / p.hbarc_MeV_fm, p; p_num=32)
    fields = Models.GasLiquidEquationSet.field_contributions(dens, p)
    masses = Models.GasLiquidEquationSet.effective_masses(st, p)
    @test fields.W ≈ p.f_omega * dens.rho_B atol=1e-12
    @test fields.R ≈ p.f_rho * dens.rho_3 atol=1e-12
    @test masses.mp ≈ p.m_nucleon_inv_fm - st.S - st.D atol=1e-12
    @test masses.mn ≈ p.m_nucleon_inv_fm - st.S + st.D atol=1e-12
    @test isfinite(dens.entropy)

    thermo = Models.GasLiquidThermodynamics.pressure_density_entropy_energy(st, 15.0 / p.hbarc_MeV_fm, p; p_num=32)
    @test thermo.entropy >= 0.0
    @test Models.GasLiquidThermodynamics.omega_identity_residual(thermo, 15.0 / p.hbarc_MeV_fm) ≈ 0.0 atol=1e-10
    @test thermo.rho_B ≈ thermo.rho_p + thermo.rho_n atol=1e-12
    @test thermo.rho_3 ≈ thermo.rho_p - thermo.rho_n atol=1e-12

    report = Models.GasLiquidThermodynamics.thermodynamic_consistency_report(st, 15.0 / p.hbarc_MeV_fm, p; p_num=32)
    @test report.finite
    @test abs(report.rho_p_error) < 1e-8
    @test abs(report.rho_n_error) < 1e-8
    @test abs(report.entropy_error) < 1e-8
    @test abs(report.identity_error) < 1e-10
end
