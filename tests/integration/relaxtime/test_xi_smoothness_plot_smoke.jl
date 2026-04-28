using Test
using JSON3

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "plot_xi_smoothness_batch.jl")
const FIXTURE_CSV = joinpath(REPO_ROOT, "tests", "fixtures", "relaxtime", "xi_smoothness", "transport_vs_xi.csv")

@testset "xi smoothness plot batch smoke" begin
    @test isfile(SCRIPT_PATH)
    @test isfile(FIXTURE_CSV)

    python_cmd = Sys.which("python") !== nothing ? "python" : (Sys.which("python3") !== nothing ? "python3" : nothing)
    has_matplotlib = python_cmd !== nothing && success(`$(python_cmd) -c "import matplotlib"`)
    if !has_matplotlib
        @test true
        return
    end

    tmp = mktempdir()
    sample_id = "itest_xi_plot_smoke"
    sample_result_dir = joinpath(tmp, "sample_results", sample_id)
    mkpath(sample_result_dir)
    sample_csv = joinpath(sample_result_dir, "transport_vs_xi.csv")
    cp(FIXTURE_CSV, sample_csv; force=true)

    manifest_path = joinpath(tmp, "run_manifest.json")
    manifest = Dict(
        "schema_version" => "v1",
        "samples" => [
            Dict(
                "sample_id" => sample_id,
                "result_csv" => sample_csv,
            ),
        ],
    )
    open(manifest_path, "w") do io
        write(io, JSON3.write(manifest))
    end

    fig_root = joinpath(tmp, "figures_root")
    fig_path = joinpath(fig_root, sample_id, "eta_over_s_vs_xi.png")
    sidecar_path = fig_path * ".provenance.json"

    run(`julia --project=. $SCRIPT_PATH --manifest $manifest_path --fig-root $fig_root`)

    @test isfile(fig_path)
    @test isfile(sidecar_path)
end
