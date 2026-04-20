using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Rotation model contract" begin
    m = Models.create_model(:Rotation)

    @test hasproperty(m.params, :a0)
    @test hasproperty(m.params, :T0_inv_fm)
    @test hasproperty(m.params, :n_cut)
    @test m.params.n_cut >= 5

    @test Models.gap_state_dim(m) == 3
    @test Models.polyakov_potential(m, 0.2, 0.3, 0.4) != 0.0

    T = 140.0 / 197.3269804
    mu = 300.0 / 197.3269804

    st = Models.solve_gap(m, T, mu; omega=0.0)
    @test st isa Models.MeanFieldState

    res = Models.RotationThermo.gap_residuals(st.phi[1], st.Phi, st.PhiBar, T, mu, 0.0, m.params)
    @test maximum(abs, res) < 1e-5

    comp = Models.omega_components(m, st, T, mu; omega=0.1)
    @test isfinite(comp.omega)
    @test isfinite(comp.therm)

    comp_core = Models.RotationThermo.omega_components(st.phi[1], st.Phi, st.PhiBar, T, mu, 0.1, m.params)
    chi_model = Models.calculate_chiral(m, st.phi; omega=0.1)
    vac_model = Models.vacuum_contribution(m, comp_core.masses; omega=0.1)
    therm_model = Models.thermal_contribution(m, comp_core.masses, st.Phi, st.PhiBar, mu, T; omega=0.1)
    @test isapprox(chi_model, comp_core.chi; rtol=1e-10, atol=1e-12)
    @test isapprox(vac_model, comp_core.vac; rtol=1e-10, atol=1e-12)
    @test isapprox(therm_model, comp_core.therm; rtol=1e-8, atol=1e-10)

    dens = Models.number_densities(m, st, T, mu; omega=0.1)
    @test all(isfinite, dens.quark)
    @test all(isfinite, dens.antiquark)
end
