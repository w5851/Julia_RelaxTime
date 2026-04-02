using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const LEGACY_SCAN_SCRIPTS = (
    "scripts/pnjl/run_tmu_scan.jl",
    "scripts/pnjl/run_dense_trho_scan.jl",
    "scripts/pnjl/run_adaptive_trho_scan.jl",
)

@testset "Wave-E compat retirement smoke" begin
    for rel_path in LEGACY_SCAN_SCRIPTS
        @test !isfile(joinpath(PROJECT_ROOT, splitpath(rel_path)...))
    end

    @test !isdefined(Main.Models, :scan_workflow_migration_map)
    @test !isdefined(Main.Models, :scan_workflow_migration_status)
end
