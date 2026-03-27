using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const AUDIT_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "config", "WorkflowConfigAudit.jl")

@testset "workflow config audit strict mode" begin
    @test isfile(AUDIT_PATH)
    include(AUDIT_PATH)
    using .WorkflowConfigAudit: build_consumption_report

    effective = Dict{String,Any}(
        "scan" => Dict{String,Any}(
            "transport" => Dict{String,Any}("resume" => true, "overwrite" => false),
            "cross_section" => Dict{String,Any}("processes" => ["ud_to_ud"]),
        ),
        "plot" => Dict{String,Any}(),
    )
    consumed = Set([
        "scan.transport.resume",
        "scan.transport.overwrite",
        "scan.cross_section.processes",
        "plot",
    ])

    rpt = build_consumption_report(effective, consumed; overridden=Set(["scan.transport.resume"]), fallback_used=false, strict=true)
    @test isempty(rpt["unused_keys"])
    @test rpt["fallback_used"] == false
    @test "scan.transport.resume" in rpt["overridden_keys"]

    consumed_incomplete = Set([
        "scan.transport.resume",
        "scan.cross_section.processes",
    ])
    @test_throws ArgumentError build_consumption_report(effective, consumed_incomplete; strict=true)
end
