using Test

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT = joinpath(ROOT, "scripts", "analysis", "pnjl_c2_limited_feasibility.jl")
const JOB_SCRIPT = joinpath(ROOT, "scripts", "analysis", "pnjl_c2_limited_feasibility_job.jl")
const WORKFLOW = joinpath(ROOT, "docs", "analysis", "governance", "diagnostic_workflow_retirement_wave2_v1", "definitions", "pnjl-c2-limited-feasibility.yml")

@testset "C2 limited feasibility contracts" begin
    script = read(SCRIPT, String)
    job = read(JOB_SCRIPT, String)
    @test occursin("pnjl_c2_limited_feasibility_v1", script)
    @test occursin("density_feasible_candidate", script)
    @test occursin("density_cap_contract_inconclusive", script)
    @test occursin("density_maxwell_candidate_inconclusive", script)
    @test occursin("density_oracle_inconclusive", script)
    @test occursin("density_performance_inconclusive", script)
    @test occursin("solver_called\" => false", script)
    @test occursin("CALCULATION_SHA_RE", script)
    @test occursin("4c9703c3be45b76608ab57d375082e29418bfd05", read(WORKFLOW, String))
    @test occursin("RHO_FINE_STEP = 0.003125", job)
    @test occursin("DENSITY_ANCHORS", job)
    @test occursin("solver_called\" => true", job)
    @test occursin("reference_write\" => false", job)
    @test occursin("fine_pool.csv", job)
    @test occursin("_field(row, :finite, false)", script)
    @test occursin("@inline function _bool", script)
    @test occursin("value isa Bool && return value", script)
end
