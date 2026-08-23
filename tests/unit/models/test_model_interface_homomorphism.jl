using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "model interface homomorphism" begin
    expected_dims = Dict(
        :NJL => 3,
        :NJL2 => 2,
        :PNJL => 5,
        :PNJLMagnetic => 5,
        :RPNJL => 5,
        :Rotation => 3,
        :GasLiquid => 4,
    )

    for (kind, expected_dim) in expected_dims
        @testset "model=$(kind)" begin
            model = kind === :PNJLMagnetic ?
                Models.create_model(kind; eB_fm2=0.1) : Models.create_model(kind)

            @test hasmethod(Models.solve_gap, Tuple{typeof(model), Any, Any})
            @test hasmethod(Models.number_densities, Tuple{typeof(model), Any, Any, Any})
            @test hasmethod(Models.calculate_mass_vec, Tuple{typeof(model), Any})
            @test hasmethod(Models.gap_state_dim, Tuple{typeof(model)})

            dim = Models.gap_state_dim(model)
            @test dim == expected_dim
        end
    end
end
