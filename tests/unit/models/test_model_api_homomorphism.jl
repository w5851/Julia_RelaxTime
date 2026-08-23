using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "models API homomorphism" begin
    for kind in (:NJL, :NJL2, :PNJL, :RPNJL, :PNJLMagnetic, :Rotation, :GasLiquid)
        model = kind === :PNJLMagnetic ?
            Models.create_model(kind; eB_fm2=0.1) : Models.create_model(kind)
        caps = Models.model_capabilities(model)
        @test hasproperty(caps, :supports_solve_gap)
        @test hasproperty(caps, :supports_model_thermo)
        @test hasproperty(caps, :supports_number_densities)
    end
end
