using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "mode accessor schema" begin
    schema = Models.schema_for_model(:PNJL)

    @testset "FixedMu 索引映射" begin
        x = [-1.84, -1.84, -2.23, 0.5, 0.5]
        sv = Models.state_view(schema, x)
        @test length(sv) == 5
        @test sv[1] == x[1]
        @test sv[5] == x[5]
    end

    @testset "FixedRho 索引映射" begin
        x = [-1.84, -1.84, -2.23, 0.5, 0.5, 0.11, 0.12, 0.13]
        sv = Models.state_view(schema, x)
        mv = Models.mu_view(schema, x; mu_dim=3)
        @test length(sv) == 5
        @test length(mv) == 3
        @test mv[1] == x[6]
        @test mv[3] == x[8]
    end

    @testset "维度错配抛 ArgumentError" begin
        @test_throws ArgumentError Models.state_view(schema, [1.0, 2.0, 3.0])
        @test_throws ArgumentError Models.mu_view(schema, [-1.0, -1.0, -2.0, 0.5, 0.5, 0.1, 0.2]; mu_dim=3)
    end
end
