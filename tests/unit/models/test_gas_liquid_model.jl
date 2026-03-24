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

    dens = Models.number_densities(m, st, T, mu)
    @test all(isfinite, dens.quark)
    @test all(isfinite, dens.antiquark)
end
