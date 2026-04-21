using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

@testset "models reverse dependency governance" begin
    tracked_files = [
        joinpath(PROJECT_ROOT, "src", "models", "derivatives", "ConservedChargeSusceptibilities.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "solver", "diff", "Targets.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "solver", "diff", "JacobianEngine.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "workflow_apps", "MesonMassWorkflow.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "precompile", "registry.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "scans", "ScanCommon.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "scans", "TmuScan.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "scans", "TrhoScan.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "workflow_apps", "TransportWorkflow.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "pnjl_physics", "core", "EquilibriumFacade.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "solver", "spec", "Conditions.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "solver", "orchestrator", "SeedStrategies.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "scans", "ScanResultFinalize.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "variants", "rotation", "workflows", "RotationWorkflow.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "variants", "gas_liquid", "workflows", "GasLiquidWorkflow.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "pnjl_physics", "PNJLMagneticModel.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "derivatives", "ThermoDerivatives.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "pnjl_physics", "core", "ModelThermodynamics.jl"),
        joinpath(PROJECT_ROOT, "src", "models", "workflow_engine", "adapters", "WorkflowAdapter.jl"),
    ]

    for file in tracked_files
        content = read(file, String)
        @test !occursin("Main.Models", content)
    end
end
