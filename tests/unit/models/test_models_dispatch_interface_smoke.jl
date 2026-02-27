using Test

_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

@testset "Models dispatch interface smoke" begin
    # Interface freeze: these generic methods must exist on AbstractQCDModel
    @test hasmethod(Models.solve_gap, Tuple{Models.AbstractQCDModel, Any, Any})
    @test hasmethod(Models.omega, Tuple{Models.AbstractQCDModel, Any, Any, Any})
    @test hasmethod(Models.omega_components, Tuple{Models.AbstractQCDModel, Any, Any, Any})
    @test hasmethod(Models.number_densities, Tuple{Models.AbstractQCDModel, Any, Any, Any})

    # Factory must continue to return AbstractQCDModel descendants.
    m_njl = Models.create_model(:NJL; profile="default")
    m_pnjl = Models.create_model(:PNJL; profile="default", physics_profile="default")
    m_rpnjl = Models.create_model(:RPNJL; profile="default", physics_profile="default")

    @test m_njl isa Models.AbstractQCDModel
    @test m_pnjl isa Models.AbstractQCDModel
    @test m_rpnjl isa Models.AbstractQCDModel
end
