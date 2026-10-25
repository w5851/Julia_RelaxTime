using Test
if !isdefined(Main,:NativeTransportRegression)
    include(joinpath(@__DIR__,"native_transport_fixture.jl"))
end

# The former v1 CSV remains historical evidence. Its underresolved sigma grid
# is no longer the production contract; retain the original numerical tolerance.
@testset "phase-guided native transport threshold regression" begin
    NativeTransportRegression.verify()
end
