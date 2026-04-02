using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "StateSchema" begin
    @testset "schema_for_model(:PNJL) contract" begin
        schema = Models.schema_for_model(:PNJL)
        @test schema.model_kind == :PNJL
        @test schema.fields == (:phi_u, :phi_d, :phi_s, :Phi, :PhiBar)
        @test Models.state_dim(schema) == 5
    end

    @testset "flatten/unflatten roundtrip" begin
        schema = Models.schema_for_model(:PNJL)
        named_state = (
            phi_u = -1.8,
            phi_d = -1.7,
            phi_s = -2.1,
            Phi = 0.25,
            PhiBar = 0.27,
        )
        x = Models.flatten_state(schema, named_state)
        @test x == [-1.8, -1.7, -2.1, 0.25, 0.27]

        state_back = Models.unflatten_state(schema, x)
        @test state_back == named_state
    end

    @testset "unflatten rejects invalid length" begin
        schema = Models.schema_for_model(:PNJL)
        @test_throws ArgumentError Models.unflatten_state(schema, [-1.0, -2.0])
    end
end
