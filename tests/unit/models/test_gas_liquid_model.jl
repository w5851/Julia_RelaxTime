using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "GasLiquid model contract" begin
    m = Models.create_model(:GasLiquid)

    @test hasproperty(m.params, :f_sigma)
    @test hasproperty(m.params, :f_omega)
    @test hasproperty(m.params, :f_rho)
    @test m.params.f_sigma ≈ 10.329 atol=1e-6

    @test Models.gap_state_dim(m) == 4

    T = 120.0 / 197.3269804
    mu = 700.0 / 197.3269804

    st = Models.solve_gap(m, T, mu)
    @test st isa Models.MeanFieldState

    eq_st = Models.GasLiquidEquationSet.GasLiquidState(st.phi[1], st.phi[2], mu, mu)
    res = Models.GasLiquidEquationSet.field_residuals(eq_st, T, m.params)
    @test maximum(abs, res) < 1e-5

    comp = Models.omega_components(m, st, T, mu)
    @test isfinite(comp.omega)
    @test comp.poly == 0.0

    gl_st = Models.GasLiquidEquationSet.GasLiquidState(st.phi[1], st.phi[2], mu, mu)
    comp_core = Models.GasLiquidThermodynamics.omega_components(gl_st, T, m.params)
    masses_model = Models.calculate_mass_vec(m, st.phi)
    chi_model = Models.calculate_chiral(m, st.phi)
    vac_model = Models.vacuum_contribution(m, masses_model; mu_vec=mu)
    therm_model = Models.thermal_contribution(m, masses_model, st.Phi, st.PhiBar, mu, T)

    @test isapprox(masses_model[1], comp_core.masses[1]; rtol=1e-10, atol=1e-12)
    @test isapprox(masses_model[2], comp_core.masses[2]; rtol=1e-10, atol=1e-12)
    @test isapprox(chi_model, comp_core.chi; rtol=1e-10, atol=1e-12)
    @test isapprox(vac_model, comp_core.vac; rtol=1e-10, atol=1e-12)
    @test isapprox(therm_model, comp_core.therm; rtol=1e-8, atol=1e-10)

    dens = Models.number_densities(m, st, T, mu)
    @test all(isfinite, dens.quark)
    @test all(isfinite, dens.antiquark)
end
