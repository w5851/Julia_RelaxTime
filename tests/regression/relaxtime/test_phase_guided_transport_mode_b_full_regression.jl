using Test
if !isdefined(Main,:NativeTransportRegression)
    include(joinpath(@__DIR__,"native_transport_fixture.jl"))
end

@testset "phase-guided native transport qualified twelve-point regression" begin
    NativeTransportRegression.verify(;full=true)
end
