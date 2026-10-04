using Test

include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "analyze_deps.jl"))

@testset "dependency audit reads current source and enforces directions" begin
    mktempdir() do root
        for dir in ("src/models/workflow_apps", "src/relaxtime", "src/simulation", "src/utils", "docs/architecture")
            mkpath(joinpath(root, dir))
        end
        write(joinpath(root, "src", "models", "Models.jl"), "")
        write(joinpath(root, "src", "relaxtime", "core.jl"), "")
        write(joinpath(root, "docs", "architecture", "dependencies.mmd"), "stale presentation, not parsed")
        write(joinpath(root, "src", "models", "workflow_apps", "Flow.jl"),
              """include(joinpath(@__DIR__, "..", "..", "relaxtime", "core.jl"))""")
        write(joinpath(root, "src", "simulation", "Server.jl"),
              """include(joinpath(@__DIR__, "..", "models", "Models.jl"))""")
        result = DependencyAudit.audit_dependencies(root)
        @test length(result.edges) == 2
        @test isempty(result.violations)
        @test isempty(result.unresolved)
        @test DependencyAudit.main(; root, strict=true, io=IOBuffer()) == 0
        write(joinpath(root, "src", "utils", "Bad.jl"),
              """include(joinpath(@__DIR__, "..", "models", "Models.jl"))""")
        @test DependencyAudit.audit_dependencies(root).violations ==
              [("src/utils/Bad.jl", "src/models/Models.jl")]
        @test DependencyAudit.main(; root, strict=true, io=IOBuffer()) == 1
        write(joinpath(root, "src", "utils", "Bad.jl"), """include("missing.jl")""")
        @test DependencyAudit.main(; root, strict=true, io=IOBuffer()) == 1
        write(joinpath(root, "src", "utils", "Bad.jl"), "include(runtime_path)")
        result = DependencyAudit.audit_dependencies(root)
        @test length(result.unresolved) == 1
        output = IOBuffer()
        DependencyAudit.write_review(output, result)
        @test occursin("include(runtime_path)", String(take!(output)))
    end
end

@testset "script bridge exceptions are exact" begin
    @test haskey(DependencyAudit.SCRIPT_BRIDGES,
        ("src/models/workflow_apps/ChargedGBUResearchWorkflow.jl",
         "scripts/analysis/relaxtime/causal_gbu_infinite_qgate.jl"))
    @test !DependencyAudit.allowed_edge("src/models/workflow_apps/New.jl", "scripts/analysis/new.jl")
    @test !DependencyAudit.allowed_edge("src/models/core/New.jl", "src/relaxtime/RelaxTime.jl")
    @test !DependencyAudit.allowed_edge("src/simulation/New.jl", "src/models/solver/Solver.jl")
end
