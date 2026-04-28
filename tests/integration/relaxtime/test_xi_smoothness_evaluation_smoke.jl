using Test
using JSON3

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "evaluate_xi_smoothness.jl")
const FIXTURE_ROOT = joinpath(REPO_ROOT, "tests", "fixtures", "relaxtime", "xi_smoothness", "eval_input")

@testset "xi smoothness evaluation smoke" begin
    @test isfile(SCRIPT_PATH)
    @test isfile(joinpath(FIXTURE_ROOT, "run_manifest.json"))

    tmp = mktempdir()
    eval_input = joinpath(tmp, "eval_input")
    cp(FIXTURE_ROOT, eval_input; force=true)

    run_manifest = joinpath(eval_input, "run_manifest.json")
    out_root = joinpath(tmp, "evaluation_out")

    run(`julia --project=. $SCRIPT_PATH --manifest $run_manifest --out-root $out_root`)

    @test isfile(joinpath(out_root, "smoothness_scores.csv"))
    @test isfile(joinpath(out_root, "smoothness_flags.csv"))
    @test isfile(joinpath(out_root, "review_queue.csv"))
end

@testset "xi smoothness evaluation skips invalid manifest samples and keeps valid ones" begin
    tmp = mktempdir()
    eval_input = joinpath(tmp, "eval_input")
    cp(FIXTURE_ROOT, eval_input; force=true)

    bad_dir = joinpath(eval_input, "sample_results", "S003")
    mkpath(bad_dir)
    bad_csv = joinpath(bad_dir, "gap_transport_scan.csv")
    write(bad_csv, join([
        "xi,tau_u,tau_s,eta_over_s,zeta_over_s,sigma_over_T",
        "-0.2,1.0,1.1,0.1,0.02,0.4",
        "0.0,bad,1.1,0.1,0.02,0.4",
        "0.2,1.0,1.1,0.1,0.02,0.4",
    ], '\n'))

    manifest_path = joinpath(eval_input, "run_manifest_mixed.json")
    manifest_payload = Dict(
        "schema_version" => "v1",
        "samples" => [
            Dict("sample_id" => "S001", "status" => "success", "result_csv" => "sample_results/S001/gap_transport_scan.csv"),
            Dict("sample_id" => "S002", "status" => "failed", "result_csv" => "sample_results/S002/gap_transport_scan.csv"),
            Dict("sample_id" => "S003", "status" => "success", "result_csv" => "sample_results/S003/gap_transport_scan.csv"),
            Dict("sample_id" => "S004", "status" => "success"),
            Dict("sample_id" => "S005", "status" => "success", "result_csv" => "sample_results/S005/gap_transport_scan.csv"),
        ],
    )
    write(manifest_path, JSON3.write(manifest_payload))

    out_root = joinpath(tmp, "evaluation_out")
    stdout_text = read(`julia --project=. $SCRIPT_PATH --manifest $manifest_path --out-root $out_root`, String)

    scores_path = joinpath(out_root, "smoothness_scores.csv")
    @test isfile(scores_path)
    @test isfile(joinpath(out_root, "smoothness_flags.csv"))
    @test isfile(joinpath(out_root, "review_queue.csv"))

    scores_text = read(scores_path, String)
    @test occursin("S001", scores_text)
    @test !occursin("S002", scores_text)
    @test !occursin("S003", scores_text)
    @test !occursin("S004", scores_text)
    @test !occursin("S005", scores_text)
    @test occursin("skipped_samples=4", stdout_text)
end
