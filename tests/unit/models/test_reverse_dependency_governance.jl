using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

@testset "models reverse dependency governance" begin
    tracked_files = [
        joinpath(PROJECT_ROOT, "src", "models", "derivatives", "ConservedChargeSusceptibilities.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "solver", "diff", "Targets.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "solver", "diff", "JacobianEngine.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "workflows", "MesonMassWorkflow.jl"),
    ]

    for file in tracked_files
        content = read(file, String)
        @test !occursin("Main.Models", content)
    end
end
