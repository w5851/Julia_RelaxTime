using Test
using TOML

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const AUTHORITY_MAP_DOC = joinpath(
    PROJECT_ROOT,
    "docs",
    "dev",
    "models-capability-map-and-dedup-table.md",
)
const AUTHORITY_MAP_CONFIG = joinpath(
    PROJECT_ROOT,
    "config",
    "governance",
    "models_authority_map.toml",
)
const AUTHORITY_MAP_CHECK_SCRIPT = joinpath(
    PROJECT_ROOT,
    "scripts",
    "dev",
    "check_models_authority_map.jl",
)

@testset "models authority map document contract" begin
    @test isfile(AUTHORITY_MAP_DOC)
    @test isfile(AUTHORITY_MAP_CONFIG)
    @test isfile(AUTHORITY_MAP_CHECK_SCRIPT)

    if isfile(AUTHORITY_MAP_DOC)
        content = read(AUTHORITY_MAP_DOC, String)
        @test occursin("## Capability Map", content)
        @test occursin("## Dedup Adjudication", content)
        @test occursin("authority boundary", content)
        @test occursin("compatibility-only", content)
    end

    if isfile(AUTHORITY_MAP_CONFIG) && isfile(AUTHORITY_MAP_DOC)
        parsed = TOML.parsefile(AUTHORITY_MAP_CONFIG)
        @test haskey(parsed, "authority_map")
        amap = parsed["authority_map"]
        @test haskey(amap, "doc")
        @test haskey(amap, "required_markers")

        doc_rel = String(amap["doc"])
        @test doc_rel == "docs/dev/models-capability-map-and-dedup-table.md"

        markers = amap["required_markers"]
        @test markers isa AbstractVector
        @test !isempty(markers)

        content = read(AUTHORITY_MAP_DOC, String)
        for marker in markers
            @test marker isa AbstractString
            @test !isempty(strip(String(marker)))
            @test occursin(String(marker), content)
        end
    end

    if isfile(AUTHORITY_MAP_CHECK_SCRIPT)
        cmd = `julia --project=$(PROJECT_ROOT) $(AUTHORITY_MAP_CHECK_SCRIPT)`
        result = run(cmd; wait=false)
        wait(result)
        @test result.exitcode == 0
    end
end
