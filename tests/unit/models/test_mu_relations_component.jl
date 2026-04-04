using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Constants_PNJL)
    include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
end

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const P = Models.pnjl_module()

@testset "mu relations component" begin
    thermal_nodes = P.cached_nodes(8, 4)
    params = P.GapParams(100.0 / 197.327, thermal_nodes, 0.0)

    @testset "equal relation in FixedRho" begin
        mode = Models.FixedRho(0.2)
        cond = P.build_conditions(mode, params)
        x = [-1.5, -1.5, -2.1, 0.2, 0.2, 1.1, 1.1, 1.1]
        r = cond([params.T_fm], x)
        @test length(r) == 8
        @test isapprox(r[6], 0.0; atol=1e-12)
        @test isapprox(r[7], 0.0; atol=1e-12)
    end

    @testset "equal relation in FixedEntropy" begin
        mode = Models.FixedEntropy(0.5)
        cond = P.build_conditions(mode, params)
        x = [-1.5, -1.5, -2.1, 0.2, 0.2, 1.2, 1.2, 1.2]
        r = cond([params.T_fm], x)
        @test length(r) == 8
        @test isapprox(r[6], 0.0; atol=1e-12)
        @test isapprox(r[7], 0.0; atol=1e-12)
    end

    @testset "asymmetric relation constraints" begin
        mode = Models.FixedAsymmetricRho(0.2, 1.0, 0.0)
        cond = P.build_conditions(mode, params)
        x = [-1.5, -1.5, -2.1, 0.2, 0.2, 1.1, 1.0, 0.9]
        r = cond([params.T_fm], x)
        @test length(r) == 8
        @test isfinite(r[6])
        @test isfinite(r[7])
        @test isfinite(r[8])
    end
end
