using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

@testset "solver include topology lock" begin
    models_path = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
    models_source = read(models_path, String)

    legacy_direct_includes = (
        "GapSolver.jl",
        "ImplicitGapLegacy.jl",
        "ConstraintSolver.jl",
        "ImplicitProblem.jl",
        "ImplicitAdapters.jl",
        "ConstraintModes.jl",
        "ConstraintComponents.jl",
        "SchemaAdapter.jl",
        "PrimaryStrategy.jl",
        "ProblemSpec.jl",
        "SolverDiagnostics.jl",
        "SolverDiagnosticsTypes.jl",
        "SolverRuntimeConfig.jl",
        "ProblemSpecOrchestrator.jl",
        "StateSchema.jl",
        "CandidateGovernance.jl",
        "SeedStrategies.jl",
        "Conditions.jl",
        "GenericRootEngine.jl",
        "WeightedFallback.jl",
        "Solver.jl",
    )

    @test occursin("include(joinpath(@__DIR__, \"solver\", \"topology.jl\"))", models_source)

    for file in legacy_direct_includes
        snippet = "include(joinpath(@__DIR__, \"solver\", \"$(file)\"))"
        @test !occursin(snippet, models_source)
    end
end
