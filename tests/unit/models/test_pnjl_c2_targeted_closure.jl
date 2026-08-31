using Test

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const JOB = joinpath(ROOT, "scripts", "analysis", "pnjl_c2_targeted_closure_job.jl")
const WORKFLOW = joinpath(ROOT, "docs", "analysis", "governance", "diagnostic_workflow_retirement_wave2_v1", "definitions", "pnjl-c2-targeted-closure-v1.yml")

@testset "C2 targeted closure workflow contract" begin
    job = read(JOB, String)
    workflow = read(WORKFLOW, String)
    @test occursin("pnjl_c2_targeted_closure_job_v1", job)
    @test occursin("3c5f6b3c9bd535cff7657364dadb2efc31f2ea48", job)
    @test occursin("target_index = findfirst", job)
    @test occursin("target = targets[target_index]", job)
    @test !occursin("target = findfirst", job)
    @test occursin("production_hybrid,independent_oracle", workflow)
    @test occursin("three_crossing_endpoint_local_v2", job)
    @test occursin("unique_three_crossing_sign_change_v2", job)
    @test occursin("oracle_labels_used_for_routing", job)
    @test occursin("reference_write", job)
    @test occursin("rerun_failed_only", workflow)
    @test occursin("source_artifact_name", workflow)
    @test occursin("cross_axis_audit", workflow)
    @test count(line -> occursin("GH_TOKEN: \${{ github.token }}", line), split(workflow, '\n')) >= 2
    @test occursin("c2-targeted-regression-", workflow)
    @test occursin("c2-targeted-cep-", workflow)
end
