using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT_PATH = joinpath(PROJECT_ROOT, "scripts", "pnjl", "run_tmu_scan.jl")

@testset "run_tmu_scan script import compatibility" begin
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
    using Main.Models: TmuScanConfig

    @test begin
        include(SCRIPT_PATH)
        true
    end
end
