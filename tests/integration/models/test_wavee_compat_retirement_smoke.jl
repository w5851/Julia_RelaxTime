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

        status = Main.Models.scan_workflow_migration_status(rel_path)
        @test status.status in (:removed, :archived)
        @test status.route == :unified_cli
        @test status.deprecation_ready == true
        @test occursin("scripts/models/run_unified_scan.jl", status.unified_entry)
    end
end
