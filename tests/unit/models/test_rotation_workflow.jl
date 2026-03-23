using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Rotation workflow entrypoints" begin
    @test isdefined(Models, :rotation_workflow_module)
    mod = Models.rotation_workflow_module()
    @test isdefined(mod, :solve_rotation_point)

    T = 150.0 / 197.3269804
    mu = 300.0 / 197.3269804
    out = Models.solve_rotation_point(T, mu; omega=0.05)

    @test out.converged
    @test isfinite(out.pressure)
    @test isfinite(out.rho)
    @test isfinite(out.entropy)
    @test isfinite(out.energy)
    @test haskey(out, :dP_domega)
end
