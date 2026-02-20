using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "core", "MagneticThermodynamics.jl"))

@testset "magnetic coupling G(B)" begin
    g0 = Main.MagneticThermodynamics.coupling_GB(0.0)
    @test isapprox(g0, Main.Constants_PNJL.G_fm2; rtol=1e-12)

    g1 = Main.MagneticThermodynamics.coupling_GB(0.05)
    g2 = Main.MagneticThermodynamics.coupling_GB(0.10)
    @test isfinite(g1)
    @test isfinite(g2)
    @test g1 > 0
    @test g2 > 0
end
