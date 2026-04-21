using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

_models_entry = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
if !(isdefined(Main, :Models) && isdefined(Main.Models, :omega) && isdefined(Main.Models, :solve_gap))
    Base.include(Main, _models_entry)
end

@testset "Models PNJL Λ injection (vacuum term)" begin
    m_default = Models.create_model(:PNJL; profile="default", physics_profile="default")
    m_bigΛ = Models.create_model(:PNJL; profile="unittest_lambda", physics_profile="default")

    @test m_default.consts.Λ_inv_fm != m_bigΛ.consts.Λ_inv_fm
    @test m_bigΛ.consts.Λ_inv_fm > m_default.consts.Λ_inv_fm

    masses = SVector{3, Float64}(0.30, 0.40, 0.50)

    vac_default = Models.vacuum_contribution(m_default, masses)
    vac_bigΛ = Models.vacuum_contribution(m_bigΛ, masses)

    @test isfinite(vac_default)
    @test isfinite(vac_bigΛ)

    # Larger cutoff => larger vacuum integral => more negative vacuum contribution.
    @test vac_bigΛ < vac_default
end
