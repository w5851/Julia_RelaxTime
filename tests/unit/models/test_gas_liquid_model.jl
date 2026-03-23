using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "GasLiquid model contract" begin
    m = Models.create_model(:GasLiquid)

    @test Models.gap_state_dim(m) == 4

    T = 120.0 / 197.3269804
    mu = 700.0 / 197.3269804

    st = Models.solve_gap(m, T, mu)
    @test st isa Models.MeanFieldState

    comp = Models.omega_components(m, st, T, mu)
    @test isfinite(comp.omega)
    @test comp.poly == 0.0

    dens = Models.number_densities(m, st, T, mu)
    @test all(isfinite, dens.quark)
    @test all(isfinite, dens.antiquark)
end
