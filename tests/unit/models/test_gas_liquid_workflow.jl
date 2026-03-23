using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "GasLiquid workflow entrypoints" begin
    @test isdefined(Models, :gas_liquid_workflow_module)
    mod = Models.gas_liquid_workflow_module()
    @test isdefined(mod, :solve_gas_liquid_point)

    T = 120.0 / 197.3269804
    mu = 700.0 / 197.3269804
    out = Models.solve_gas_liquid_point(T, mu)

    @test haskey(out, :pressure)
    @test haskey(out, :rho)
    @test haskey(out, :entropy)
    @test haskey(out, :energy)
    @test out.converged == true
    @test isfinite(out.pressure)
    @test isfinite(out.rho)
end
