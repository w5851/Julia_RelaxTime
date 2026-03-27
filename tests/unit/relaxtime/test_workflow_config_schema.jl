using Test
using TOML

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

const SCHEMA_PATH = joinpath(REPO_ROOT, "config", "workflows", "relaxtime", "schema", "relaxtime_workflow_schema_v1.toml")
const ALIASES_PATH = joinpath(REPO_ROOT, "config", "workflows", "relaxtime", "schema", "aliases_v1.toml")
const DEFAULT_PATH = joinpath(REPO_ROOT, "config", "workflows", "relaxtime", "default.toml")
const PROFILE_TRANSPORT_PATH = joinpath(REPO_ROOT, "config", "workflows", "relaxtime", "profiles", "muB0_transport_xi.toml")
const PROFILE_XS_PATH = joinpath(REPO_ROOT, "config", "workflows", "relaxtime", "profiles", "muB0_cross_section_xi.toml")

@testset "relaxtime workflow schema assets" begin
    @test isfile(SCHEMA_PATH)
    @test isfile(ALIASES_PATH)
    @test isfile(DEFAULT_PATH)
    @test isfile(PROFILE_TRANSPORT_PATH)
    @test isfile(PROFILE_XS_PATH)
end

@testset "relaxtime workflow defaults and frozen contracts" begin
    cfg = TOML.parsefile(DEFAULT_PATH)

    @test get(cfg, "schema_version", "") == "v1"

    scan = get(cfg, "scan", Dict{String, Any}())
    transport = get(scan, "transport", Dict{String, Any}())
    xs = get(scan, "cross_section", Dict{String, Any}())

    @test get(transport, "muB_MeV", nothing) == 0.0
    @test get(transport, "xi_list", Any[]) == Any[-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]

    @test get(xs, "muB_MeV", nothing) == 0.0
    @test get(xs, "T_list_MeV", Any[]) == Any[150.0, 250.0]
    @test get(xs, "processes", Any[]) == Any["ud_to_ud", "us_to_us", "udbar_to_udbar", "usbar_to_usbar"]

    energy = get(xs, "energy", Dict{String, Any}())
    @test haskey(energy, "mode")
end

@testset "relaxtime aliases include process compatibility" begin
    alias_cfg = TOML.parsefile(ALIASES_PATH)
    proc_aliases = get(get(alias_cfg, "process_aliases", Dict{String, Any}()), "mappings", Dict{String, Any}())

    @test get(proc_aliases, "udtoud", "") == "ud_to_ud"
    @test get(proc_aliases, "ustous", "") == "us_to_us"
    @test get(proc_aliases, "udbartoudbar", "") == "udbar_to_udbar"
    @test get(proc_aliases, "usbartousbar", "") == "usbar_to_usbar"
end
