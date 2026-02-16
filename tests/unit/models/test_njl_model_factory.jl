using Test

const _MODELS_ENTRY = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_MODELS_ENTRY)
end
using .Models

using StaticArrays

@testset "Models factory: NJL" begin
    m = Models.create_model(:NJL; profile="default")
    @test m isa Models.NJLModel
    @test m.params.label == "njl-standard"

    φ = @SVector [0.01, 0.02, 0.03]
    masses = Models.calculate_mass_vec(m, φ)
    @test masses == Models.NJLCore.calculate_mass_vec(m.params, φ)

    χ = Models.calculate_chiral(m, φ)
    @test χ == Models.NJLCore.chiral_potential(m.params, φ)

    @test Models.polyakov_potential(m, 0.1, 0.2, 0.3) == 0.0
end
