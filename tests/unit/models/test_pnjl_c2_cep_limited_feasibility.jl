using Test

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const EVALUATOR = joinpath(ROOT, "scripts", "analysis", "pnjl_c2_cep_limited_feasibility.jl")
const JOB = joinpath(ROOT, "scripts", "analysis", "pnjl_c2_cep_limited_feasibility_job.jl")
const WORKFLOW = joinpath(ROOT, ".github", "workflows", "pnjl-c2-cep-limited-feasibility.yml")

@testset "C2 CEP limited feasibility contracts" begin
    evaluator = read(EVALUATOR, String)
    job = read(JOB, String)
    workflow = read(WORKFLOW, String)
    @test occursin("pnjl_c2_cep_limited_feasibility_v2", evaluator)
    @test occursin("pnjl_c2_cep_limited_feasibility_job_v2", evaluator)
    @test occursin("_production_maxwell_options", evaluator)
    @test occursin("_production_maxwell_or_empty", evaluator)
    @test occursin("hybrid_states=endpoints.hybrid_states", evaluator)
    @test occursin("point_request_reconciliation", evaluator)
    @test occursin("oracle_labels_used_for_routing=false", job)
    @test occursin("point_requests=item.cache.unique_solves + item.cache.cache_hits", job) ||
        occursin("point_requests == unique_solves + cache_hits", job)
    @test occursin("method=\"hybrid\"", job)
    @test occursin("method=\"oracle\"", job)
    @test occursin("pnjl_c2_cep_limited_feasibility_v2", workflow)
    @test occursin("options: [cep]", workflow)
    @test occursin("max-parallel: 17", workflow)
    @test occursin("rerun_failed_only", workflow)
end
