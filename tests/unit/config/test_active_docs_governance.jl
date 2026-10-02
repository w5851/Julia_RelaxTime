using Dates
using Test

if !isdefined(Main, :ActiveDocsGovernance)
    include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "check_active_docs_governance.jl"))
end

@testset "Active document review follows status, not age" begin
    mktempdir() do root
        active = joinpath(root, "docs", "dev", "active")
        mkpath(active)
        write(joinpath(active, "2026-01-01_long_running_task.md"), "# Task\nStatus: active\n")
        result = ActiveDocsGovernance.review_documents(root; current_date=Date(2026, 10, 2))
        @test isempty(result.violations)
        @test length(result.advisories) == 1
        @test isfile(joinpath(active, "2026-01-01_long_running_task.md"))

        write(joinpath(active, "2026-09-30_recent_task.md"), "# Task\n")
        @test length(ActiveDocsGovernance.review_documents(root; current_date=Date(2026, 10, 2)).advisories) == 1
        write(joinpath(active, "2026-99-99_invalid_date.md"), "# Task\n")
        @test length(ActiveDocsGovernance.review_documents(root; current_date=Date(2026, 10, 2)).violations) == 1
    end
end
