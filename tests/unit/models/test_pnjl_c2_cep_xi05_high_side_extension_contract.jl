using Test

const HIGH_SIDE_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const HIGH_SIDE_JOB = joinpath(HIGH_SIDE_ROOT, "scripts", "analysis",
    "pnjl_c2_cep_xi05_high_side_extension_job.jl")
const HIGH_SIDE_PLAN = joinpath(HIGH_SIDE_ROOT, "docs", "analysis", "pnjl",
    "c2_followups", "c2_cep_xi05_high_side_extension_v1", "temperature_plan.csv")

@testset "C2 xi=0.5 high-side extension contract" begin
    source = read(HIGH_SIDE_JOB, String)
    @test occursin("HIGH_SIDE_CALCULATION_SHA", source)
    @test occursin("HIGH_SIDE_ANCHOR_T = 107.0625", source)
    @test occursin("HIGH_SIDE_STEP = 0.0625", source)
    @test occursin("HIGH_SIDE_TEMPERATURES = (107.125, 107.1875, 107.25)", source)
    @test occursin("oracle_labels_used_for_routing=false", source)
    @test occursin("_validate_high_cost_rows(cost_rows)", source)
    @test occursin("point_requests == item.unique_solves + item.cache_hits", source)
    @test !occursin("oracle_status", source)

    rows = readlines(HIGH_SIDE_PLAN)
    @test length(rows) == 4
    @test rows[2] == "0.5,1,high_extension_1,107.125,107.0625,0.0625,high"
    @test rows[3] == "0.5,2,high_extension_2,107.1875,107.0625,0.0625,high"
    @test rows[4] == "0.5,3,high_extension_3,107.25,107.0625,0.0625,high"
end
