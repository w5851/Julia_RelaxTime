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

    dens = Models.number_densities(m, st, T, mu; omega=0.1)
    @test all(isfinite, dens.quark)
    @test all(isfinite, dens.antiquark)
end
