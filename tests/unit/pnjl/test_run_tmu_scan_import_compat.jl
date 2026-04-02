using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT_PATH = joinpath(PROJECT_ROOT, "scripts", "pnjl", "run_tmu_scan.jl")

@testset "run_tmu_scan script import compatibility" begin
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
    using Main.Models: TmuScanConfig

    if !isfile(SCRIPT_PATH)
        @test !isdefined(Main.Models, :scan_workflow_migration_status)
        return
    end

    @test include(SCRIPT_PATH) === nothing
end
