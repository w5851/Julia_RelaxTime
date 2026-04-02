using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "model skeleton homomorphism smoke" begin
    @testset "all model kinds expose dedicated workflow module" begin
        kinds = (:NJL, :NJL2, :PNJL, :PNJLMagnetic, :RPNJL, :Rotation, :GasLiquid)
        for kind in kinds
            mod = Main.Models.workflow_module_for(kind)
            @test mod isa Module
            if kind === :Rotation
                @test isdefined(mod, :solve_rotation_point)
            elseif kind === :GasLiquid
                @test isdefined(mod, :solve_gas_liquid_point)
            else
                @test isdefined(mod, :solve_gap_and_transport)
            end
        end
    end
end
