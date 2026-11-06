using Test

const SOP_CHECKER_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "check_sop_governance.jl"))

if !isdefined(Main, :SopGovernance)
    include(SOP_CHECKER_PATH)
end

function _write_text(path::AbstractString, content::AbstractString)
    mkpath(dirname(path))
    write(path, content)
    return path
end

_sop_content() = "# Fixture workflow\n\nInstructions for this workflow.\n"

function _registry_text(; last_verified::String="2026-07-10", duplicate_authority::Bool=false)
    second = duplicate_authority ? """

[[sop]]
id = "duplicate"
path = "docs/guides/sop/duplicate.md"
status = "active"
authoritative_for = ["fixture_workflow"]
entrypoints = ["scripts/run_fixture.jl"]
stable_entrypoints = ["scripts/run_fixture.jl"]
configs = []
last_verified = "$(last_verified)"
verification_commands = ["julia --project=. scripts/run_fixture.jl --help"]
""" : ""

    return """
schema_version = "v1"
allowed_statuses = ["draft", "active", "deprecated"]
stable_entrypoint_index = "docs/guides/scripts/README.md"
forbidden_patterns = ["src/pnjl/"]

[[sop]]
id = "fixture"
path = "docs/guides/sop/fixture.md"
status = "active"
authoritative_for = ["fixture_workflow"]
entrypoints = ["scripts/run_fixture.jl"]
stable_entrypoints = ["scripts/run_fixture.jl"]
configs = []
last_verified = "$(last_verified)"
verification_commands = ["julia --project=. scripts/run_fixture.jl --help"]
$(second)
"""
end

function _make_fixture(root::AbstractString; last_verified="2026-07-10", duplicate_authority=false)
    _write_text(joinpath(root, "scripts", "run_fixture.jl"), "println(\"fixture\")\n")
    _write_text(
        joinpath(root, "docs", "guides", "scripts", "README.md"),
        "# Fixture scripts\n\n[上级说明](../README.md)\n\nscripts\\run_fixture.jl\n",
    )
    _write_text(
        joinpath(root, "docs", "guides", "sop", "fixture.md"),
        _sop_content(),
    )
    if duplicate_authority
        _write_text(
            joinpath(root, "docs", "guides", "sop", "duplicate.md"),
            _sop_content(),
        )
    end
    _write_text(
        joinpath(root, "config", "governance", "docs_authority_map.toml"),
        _registry_text(; last_verified=last_verified, duplicate_authority=duplicate_authority),
    )
    return root
end

@testset "SOP governance registry" begin
    @test isfile(SOP_CHECKER_PATH)
    @test Main.SopGovernance.normalize_text_slashes("docs\\guides\n[上级](../README.md)") ==
        "docs/guides\n[上级](../README.md)"

    mktempdir() do temp_root
        valid_root = _make_fixture(joinpath(temp_root, "valid"))
        @test isempty(Main.SopGovernance.validate_registry(
            valid_root,
        ))

        invalid_date = _make_fixture(joinpath(temp_root, "invalid-date"); last_verified="not-a-date")
        @test any(item -> occursin("last_verified", item), Main.SopGovernance.validate_registry(invalid_date))

        no_date_root = _make_fixture(joinpath(temp_root, "no-date"))
        registry = joinpath(no_date_root, "config", "governance", "docs_authority_map.toml")
        write(registry, replace(read(registry, String), r"last_verified = [^\n]*\n" => ""))
        @test isempty(Main.SopGovernance.validate_registry(no_date_root))

        duplicate_root = _make_fixture(
            joinpath(temp_root, "duplicate");
            duplicate_authority=true,
        )
        duplicate = Main.SopGovernance.validate_registry(
            duplicate_root,
        )
        @test any(item -> occursin("is claimed by both", item), duplicate)

        rm(joinpath(valid_root, "scripts", "run_fixture.jl"))
        @test any(item -> occursin("does not exist", item), Main.SopGovernance.validate_registry(
            valid_root,
        ))
    end
end
