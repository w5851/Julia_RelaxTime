using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Solver schema adapter" begin
    @testset "registry lookup" begin
        reg = Models.SchemaRegistry()
        schema = Models.VarSchema((:x1, :x2), (:T, :mu))
        Models.register_schema!(reg, :PNJL, :fixed_mu, schema)

        @test Models.schema_for(reg, :PNJL, :fixed_mu) == schema
        @test_throws ArgumentError Models.schema_for(reg, :PNJL, :fixed_rho)
    end

    @testset "named vector roundtrip" begin
        schema = Models.VarSchema((:phi_u, :phi_d, :Phi), (:T, :mu))

        x_named = (phi_u=-1.8, phi_d=-1.7, Phi=0.2)
        theta_named = (T=0.8, mu=0.1)

        x_vec = Models.named_to_vec(x_named, schema, :x)
        theta_vec = Models.named_to_vec(theta_named, schema, :theta)

        @test x_vec == [-1.8, -1.7, 0.2]
        @test theta_vec == [0.8, 0.1]

        @test Models.vec_to_named(x_vec, schema, :x) == x_named
        @test Models.vec_to_named(theta_vec, schema, :theta) == theta_named
    end

    @testset "schema validation" begin
        schema = Models.VarSchema((:x1, :x2, :x3), (:t1, :t2))
        @test Models.validate_schema(schema; x_dim=3, theta_dim=2) === nothing
        @test_throws ArgumentError Models.validate_schema(schema; x_dim=2)
        @test_throws ArgumentError Models.validate_schema(schema; theta_dim=1)

        bad_schema = Models.VarSchema((:x1, :x1), (:t1, :t2))
        @test_throws ArgumentError Models.validate_schema(bad_schema)
    end
end
