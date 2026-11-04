using Test
if !isdefined(Main,:TransportValidationSelection)
    include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","TransportValidationSelection.jl"))
end
using .TransportValidationSelection: requested_validation_rows

@testset "Explicit transport validation point selection" begin
    source=[(point_id="T120_muB900_xi-0.5",T_MeV=120.0),
        (point_id="T160_muB450_xi0.0",T_MeV=160.0),
        (point_id="T200_muB0_xi0.5",T_MeV=200.0)]
    ids=[source[3].point_id,source[1].point_id]
    @test requested_validation_rows(source,ids)==source[[1,3]]
    @test requested_validation_rows(source,[source[2].point_id])==source[2:2]
    @test_throws ArgumentError requested_validation_rows(source,String[])
    @test_throws ArgumentError requested_validation_rows(source,[source[1].point_id,source[1].point_id])
    @test_throws ArgumentError requested_validation_rows(source,["missing"])
    @test_throws ArgumentError requested_validation_rows(vcat(source,source[1:1]),ids)
    @test_throws ArgumentError requested_validation_rows(source,[source[1].point_id,"missing"])
end
