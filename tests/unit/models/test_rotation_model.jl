using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Rotation model contract" begin
    m = Models.create_model(:Rotation)

    @test Models.gap_state_dim(m) == 3
    @test Models.polyakov_potential(m, 0.2, 0.3, 0.4) == 0.0

    T = 140.0 / 197.3269804
    mu = 300.0 / 197.3269804

    st = Models.solve_gap(m, T, mu; omega=0.0)
    @test st isa Models.MeanFieldState

    comp = Models.omega_components(m, st, T, mu; omega=0.1)
    @test isfinite(comp.omega)
    @test isfinite(comp.therm)

    dens = Models.number_densities(m, st, T, mu; omega=0.1)
    @test all(isfinite, dens.quark)
    @test all(isfinite, dens.antiquark)
end
