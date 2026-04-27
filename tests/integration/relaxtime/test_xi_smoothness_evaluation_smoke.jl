using Test

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
