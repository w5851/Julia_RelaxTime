using Test
include(normpath(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","compute_freezeout_trho_overlay.jl")))
const FTO=FreezeoutTRhoOverlay
@testset "Freezeout T-rho projection: net charges and historical domain" begin
    c=FTO.coordinates([0.224,0.256,0.],0.16)
    @test c.rho_norm≈1.0
    @test c.rho_B≈0.16
    @test c.rho_Q/c.rho_B≈0.4
    @test c.rho_S==0
    @test FTO.coordinates([-0.224,-0.256,0.],0.16).rho_norm≈-1.0
    @test FTO.coverage(120.,0.05)=="inside_historical_grid"
    @test FTO.coverage(220.,1.)=="inside_historical_grid"
    @test FTO.coverage(80.,0.1)=="temperature_outside"
    @test FTO.coverage(160.,0.01)=="density_outside"
    @test FTO.coverage(80.,0.01)=="temperature_outside;density_outside"
    @test FTO.coverage(NaN,0.1)=="nonfinite"
    @test_throws ArgumentError FTO.coordinates([1.,2.],0.16)
    @test_throws ArgumentError FTO.coordinates([1.,2.,NaN],0.16)
    @test_throws ArgumentError FTO.coordinates([1.,2.,3.],0.)
end
