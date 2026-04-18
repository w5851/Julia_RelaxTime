using Test
using JSON3

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const WORKFLOW_PATH = joinpath(REPO_ROOT, "scripts", "relaxtime", "run_manual_relaxation_scan_workflow.jl")

function _first_csv_header(path::String)
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            return [strip(x) for x in split(s, ',')]
        end
    end
    return String[]
end

@testset "manual relaxation workflow provenance smoke" begin
    @test isfile(WORKFLOW_PATH)

    base_output_dir = mktempdir()
    cmd = `julia --project=. $WORKFLOW_PATH --sections plan_a,plan_b --no-plots --overwrite --base-output-dir $base_output_dir --plan-a-Tmin 150 --plan-a-Tmax 150 --plan-a-Tstep 10 --plan-b-T-list 190 --xi-min -0.5 --xi-max -0.5 --xi-step 0.1 --mode finite_15 --tau-p-nodes 8 --tau-angle-nodes 2 --tau-phi-nodes 2 --tau-n-sigma 4 --sigma-grid-n 32 --no-bulk`
    run(cmd)

    plan_a_dir = joinpath(base_output_dir, "results", "relaxtime", "plan_a")
    plan_b_dir = joinpath(base_output_dir, "results", "relaxtime", "plan_b")

    plan_a_csv = joinpath(plan_a_dir, "gap_transport_vs_T_muB0_xi0.csv")
    plan_b_csv = joinpath(plan_b_dir, "transport_vs_xi_T190_muB0.csv")
    plan_b_merged = joinpath(plan_b_dir, "plan_b_merged.csv")

    @test isfile(plan_a_csv)
    @test isfile(plan_b_csv)
    @test isfile(plan_b_merged)

    @test isfile(joinpath(plan_a_dir, "effective_config.json"))
    @test isfile(joinpath(plan_a_dir, "run_manifest.json"))
    @test isfile(joinpath(plan_b_dir, "effective_config.json"))
    @test isfile(joinpath(plan_b_dir, "run_manifest.json"))

    merged_text = read(plan_b_merged, String)
    @test occursin("run_id", merged_text)
    @test occursin("source_file", merged_text)
    @test occursin("source_T_MeV", merged_text)

    manifest_text = read(joinpath(plan_b_dir, "run_manifest.json"), String)
    @test occursin("\"artifacts\"", manifest_text)
    @test occursin("\"summary\"", manifest_text)
    @test occursin("\"fallback_count\"", manifest_text)

    manifest_obj = JSON3.read(manifest_text)
    cwd_path = String(manifest_obj["cwd"])
    project_path = String(manifest_obj["project_path"])
    @test !occursin('\\', cwd_path)
    @test !occursin('\\', project_path)
    @test !occursin(r"^[A-Za-z]:", cwd_path)
    @test project_path == "."

    artifacts = manifest_obj["artifacts"]
    @test length(artifacts) > 0
    first_artifact_path = String(artifacts[1]["path"])
    @test !occursin('\\', first_artifact_path)
    @test !occursin(r"^[A-Za-z]:", first_artifact_path)
end
