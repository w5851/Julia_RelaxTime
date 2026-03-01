using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()
const PNJL = Models.pnjl_module()
const coupling_GB = getproperty(PNJL, :coupling_GB)
const Constants_PNJL = getproperty(PNJL, :Constants_PNJL)
const G_fm2 = getproperty(Constants_PNJL, :G_fm2)

@testset "magnetic coupling G(B)" begin
    g0 = coupling_GB(0.0)
    @test isapprox(g0, G_fm2; rtol=1e-12)

    g1 = coupling_GB(0.05)
    g2 = coupling_GB(0.10)
    @test isfinite(g1)
    @test isfinite(g2)
    @test g1 > 0
    @test g2 > 0
end
