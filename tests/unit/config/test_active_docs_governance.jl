using Dates
using Test

if !isdefined(Main, :ActiveDocsGovernance)
    include(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "check_active_docs_governance.jl"))
end

@testset "A clean checkout may have no active documents" begin
    mktempdir() do root
        mkpath(joinpath(root, "docs", "dev"))
        active = joinpath(root, "docs", "dev", "active")
        result = ActiveDocsGovernance.review_documents(root; current_date=Date(2026, 11, 8))
        @test isempty(result.violations)
        @test isempty(result.advisories)
        @test !ispath(active)
        mkpath(active)
        @test ActiveDocsGovernance.review_documents(root; current_date=Date(2026, 11, 8)) == result
    end
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

@testset "Retired directory gate blocks new files" begin
    mktempdir() do root
        notes = joinpath(root, "docs", "notes")
        mkpath(joinpath(notes, "pnjl"))
        write(joinpath(notes, "README.md"), "# notes\n")
        write(joinpath(notes, "pnjl", "kept.md"), "# kept\n")
        allowlist = joinpath(root, "retired_docs_allowlist.toml")
        write(allowlist, """
        schema_version = 1

        [directories."docs/notes"]
        status = "retired_readonly"
        allowed_files = ["README.md", "pnjl/kept.md"]
        """)
        clean = ActiveDocsGovernance.review_retired_directories(root; allowlist_path=allowlist)
        @test isempty(clean.violations)
        @test isempty(clean.advisories)

        write(joinpath(notes, "new_note.md"), "# new\n")
        blocked = ActiveDocsGovernance.review_retired_directories(root; allowlist_path=allowlist)
        @test length(blocked.violations) == 1
        @test occursin("docs/notes/new_note.md", only(blocked.violations))

        rm(joinpath(notes, "pnjl", "kept.md"))
        missing = ActiveDocsGovernance.review_retired_directories(root; allowlist_path=allowlist)
        @test any(item -> occursin("kept.md", item), missing.advisories)
    end
end

@testset "Retired gate stays off outside the repository root" begin
    mktempdir() do root
        result = ActiveDocsGovernance.review_documents(root; current_date=Date(2026, 11, 8))
        @test isempty(result.violations)
        @test isempty(result.advisories)
    end
end

@testset "Broken or missing allowlist fails closed" begin
    mktempdir() do root
        notes = joinpath(root, "docs", "notes")
        mkpath(notes)
        write(joinpath(notes, "README.md"), "# notes\n")

        broken = joinpath(root, "broken_allowlist.toml")
        write(broken, "schema_version = 1\n[directories.\"docs/notes\"\n")
        result = ActiveDocsGovernance.review_retired_directories(root; allowlist_path=broken)
        @test length(result.violations) == 1
        @test occursin("invalid retired docs allowlist", only(result.violations))

        unversioned = joinpath(root, "unversioned_allowlist.toml")
        write(unversioned, """
        schema_version = 99

        [directories."docs/notes"]
        status = "retired_readonly"
        allowed_files = ["README.md"]
        """)
        result = ActiveDocsGovernance.review_retired_directories(root; allowlist_path=unversioned)
        @test length(result.violations) == 1
        @test occursin("schema_version", only(result.violations))

        no_status = joinpath(root, "no_status_allowlist.toml")
        write(no_status, """
        schema_version = 1

        [directories."docs/notes"]
        allowed_files = ["README.md"]
        """)
        result = ActiveDocsGovernance.review_retired_directories(root; allowlist_path=no_status)
        @test length(result.violations) == 1
        @test occursin("missing status", only(result.violations))
    end
end
