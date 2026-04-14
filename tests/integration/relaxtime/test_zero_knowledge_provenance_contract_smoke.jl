using Test

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

    manifest_text = read(manifest_path, String)
    @test occursin("\"run_id\"", manifest_text)
    @test occursin("\"artifacts\"", manifest_text)
    @test occursin("\"sha256\"", manifest_text)
    @test occursin("\"rows\"", manifest_text)
end
