using Test
using JSON3

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BATCH_SCRIPT = joinpath(REPO_ROOT, "scripts", "relaxtime", "run_xi_smoothness_batch.jl")

function _write_params_csv(path::String, rows::Vector{String})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "sample_id,source,T_MeV,muq_MeV,muB_MeV,anchor_type,anchor_T_MeV,anchor_muq_MeV,delta_T,delta_muq,rng_seed")
        for row in rows
            println(io, row)
        end
    end
end

function _read_manifest(out_root::String)
    manifest_path = joinpath(out_root, "run_manifest.json")
    @test isfile(manifest_path)
    return JSON3.read(read(manifest_path, String))
end

@testset "xi smoothness batch runner dry-run manifest smoke" begin
    @test isfile(BATCH_SCRIPT)

    tmp = mktempdir()
    params_csv = joinpath(tmp, "params.csv")
    out_root = joinpath(tmp, "batch_out")
    _write_params_csv(params_csv, [
        "S001,random_uniform,150.0,50.0,150.0,none,,,,,42",
        "S002,near_phase_line,190.0,120.0,360.0,boundary,188.0,118.0,2.0,2.0,42",
    ])

    run(`julia --project=. $BATCH_SCRIPT --params $params_csv --out-root $out_root --dry-run`)
    manifest = _read_manifest(out_root)
    samples = manifest["samples"]
    @test length(samples) == 2
    for sample in samples
        @test haskey(sample, "status")
        @test haskey(sample, "result_csv")
        @test haskey(sample, "failed_points_path")
    end
end

@testset "xi smoothness batch runner dry-run overwrite keeps existing artifacts" begin
    tmp = mktempdir()
    params_csv = joinpath(tmp, "params.csv")
    out_root = joinpath(tmp, "batch_out")
    _write_params_csv(params_csv, [
        "S001,random_uniform,150.0,50.0,150.0,none,,,,,42",
    ])

    sample_dir = joinpath(out_root, "S001")
    mkpath(sample_dir)
    result_csv = joinpath(sample_dir, "gap_transport_scan.csv")
    failed_csv = joinpath(sample_dir, "failed_points.csv")
    write(result_csv, "existing_result")
    write(failed_csv, "existing_failed")

    run(`julia --project=. $BATCH_SCRIPT --params $params_csv --out-root $out_root --dry-run --overwrite`)

    @test isfile(result_csv)
    @test isfile(failed_csv)
    @test read(result_csv, String) == "existing_result"
    @test read(failed_csv, String) == "existing_failed"
end

@testset "xi smoothness batch runner dry-run takes precedence over resume" begin
    tmp = mktempdir()
    params_csv = joinpath(tmp, "params.csv")
    out_root = joinpath(tmp, "batch_out")
    _write_params_csv(params_csv, [
        "S001,random_uniform,150.0,50.0,150.0,none,,,,,42",
    ])

    sample_dir = joinpath(out_root, "S001")
    mkpath(sample_dir)
    result_csv = joinpath(sample_dir, "gap_transport_scan.csv")
    write(result_csv, "already_exists")

    run(`julia --project=. $BATCH_SCRIPT --params $params_csv --out-root $out_root --dry-run --resume`)
    manifest = _read_manifest(out_root)
    sample = manifest["samples"][1]
    @test String(sample["status"]) == "skipped"
    @test String(sample["reason"]) == "dry_run"
end

@testset "xi smoothness batch runner marks resume_hit as success when result exists" begin
    tmp = mktempdir()
    params_csv = joinpath(tmp, "params.csv")
    out_root = joinpath(tmp, "batch_out")
    _write_params_csv(params_csv, [
        "S001,random_uniform,150.0,50.0,150.0,none,,,,,42",
    ])

    sample_dir = joinpath(out_root, "S001")
    mkpath(sample_dir)
    result_csv = joinpath(sample_dir, "gap_transport_scan.csv")
    write(result_csv, "xi,tau_u,tau_s,eta_over_s,zeta_over_s,sigma_over_T\n0.0,1,1,1,1,1\n")

    run(`julia --project=. $BATCH_SCRIPT --params $params_csv --out-root $out_root --resume`)
    manifest = _read_manifest(out_root)
    sample = manifest["samples"][1]
    @test String(sample["status"]) == "success"
    @test String(sample["reason"]) == "resume_hit"
end

@testset "xi smoothness batch runner rejects duplicate sample_id" begin
    tmp = mktempdir()
    params_csv = joinpath(tmp, "params.csv")
    out_root = joinpath(tmp, "batch_out")
    _write_params_csv(params_csv, [
        "S001,random_uniform,150.0,50.0,150.0,none,,,,,42",
        "S001,near_phase_line,190.0,120.0,360.0,boundary,188.0,118.0,2.0,2.0,42",
    ])

    err = IOBuffer()
    proc = run(pipeline(
        ignorestatus(`julia --project=. $BATCH_SCRIPT --params $params_csv --out-root $out_root --dry-run`);
        stdout=devnull,
        stderr=err,
    ))
    @test proc.exitcode != 0

    msg = String(take!(err))
    @test occursin("duplicate sample_id", msg)
    @test occursin("S001", msg)
    @test occursin("line", msg)
end
