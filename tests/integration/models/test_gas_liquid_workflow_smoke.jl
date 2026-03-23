using Test

_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

@testset "GasLiquid workflow smoke" begin
    T_points = (100.0 / 197.3269804, 120.0 / 197.3269804)
    mu = 700.0 / 197.3269804

    for T in T_points
        out = Models.solve_gas_liquid_point(T, mu)
        @test out.converged
        @test isfinite(out.pressure)
        @test isfinite(out.rho)
        @test isfinite(out.entropy)
        @test isfinite(out.energy)
    end
end
