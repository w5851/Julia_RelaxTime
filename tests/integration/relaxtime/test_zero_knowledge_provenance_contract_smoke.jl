using Test
using JSON3

const _PROVENANCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "provenance_metadata.jl"))
if !isdefined(Main, :ProvenanceMetadata)
    Base.include(Main, _PROVENANCE_PATH)
end

using Main.ProvenanceMetadata

@testset "provenance sidecars minimal contract" begin
    tmp = mktempdir()
    outdir = joinpath(tmp, "out")
    mkpath(outdir)
    csv_path = joinpath(outdir, "demo.csv")

    open(csv_path, "w") do io
        println(io, "# schema: scan_csv_v1")
        println(io, "A,converged")
        println(io, "1,true")
        println(io, "2,false")
    end

    stats = ProvenanceMetadata.csv_stats(csv_path; converged_col="converged")
    @test stats.points_total == 2
    @test stats.success_count == 1
    @test stats.error_count == 1

    ctx = ProvenanceMetadata.new_run_context("scripts/relaxtime/run_gap_transport_scan.jl", ["--tmin", "120"])
    @test !isempty(ctx.run_id)
    @test !isempty(ctx.timestamp_utc)

    ProvenanceMetadata.write_run_sidecars(
        outdir;
        ctx=ctx,
        effective_config=Dict{String,Any}("demo" => true),
        artifacts=[csv_path],
        summary=Dict{String,Any}("points_total" => 2),
    )

    effective_path = joinpath(outdir, "effective_config.json")
    manifest_path = joinpath(outdir, "run_manifest.json")
    @test isfile(effective_path)
    @test isfile(manifest_path)

    manifest_obj = JSON3.read(read(manifest_path, String))
    effective_obj = JSON3.read(read(effective_path, String))
    @test String(manifest_obj["run_id"]) == ctx.run_id
    @test haskey(manifest_obj, "artifacts")
    @test length(manifest_obj["artifacts"]) == 1
    @test haskey(manifest_obj["artifacts"][1], "sha256")
    @test haskey(manifest_obj["artifacts"][1], "rows")
    @test String(manifest_obj["project_path"]) == "."
    @test !isabspath(String(manifest_obj["cwd"]))
    @test String(effective_obj["run_id"]) == ctx.run_id

    @test_throws ArgumentError ProvenanceMetadata.write_run_sidecars(
        outdir;
        ctx=ctx,
        effective_config=Dict{String,Any}("bad" => NaN),
        artifacts=[csv_path],
        summary=Dict{String,Any}("points_total" => 2),
    )
end
