using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "models dimension-agnostic solver regression" begin
    expected_dims = Dict(
        :NJL => 3,
        :NJL2 => 2,
        :PNJL => 5,
        :PNJLMagnetic => 5,
        :RPNJL => 5,
        :Rotation => 3,
        :GasLiquid => 4,
    )

    for (kind, dim) in expected_dims
        model = kind === :PNJLMagnetic ?
            Models.create_model(kind; eB_fm2=0.1) : Models.create_model(kind)
        @test Models.gap_state_dim(model) == dim

        schema = Models.schema_for_model(kind)
        @test schema.model_kind == kind
        @test Models.state_dim(schema) == dim

        named_state = NamedTuple{schema.fields}(ntuple(i -> -0.1 * i, dim))
        x = Models.flatten_state(schema, named_state)
        @test length(x) == dim
        @test Models.unflatten_state(schema, x) == named_state
    end
end
