using Test

_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

@testset "Rotation workflow smoke" begin
    mu = 300.0 / 197.3269804
    omega = 0.05
    T_points = (140.0 / 197.3269804, 160.0 / 197.3269804)

    for T in T_points
        out = Models.solve_rotation_point(T, mu; omega=omega)
        @test out.converged
        @test isfinite(out.pressure)
        @test isfinite(out.rho)
        @test isfinite(out.entropy)
        @test isfinite(out.energy)
        @test isfinite(out.dP_domega)
    end
end
