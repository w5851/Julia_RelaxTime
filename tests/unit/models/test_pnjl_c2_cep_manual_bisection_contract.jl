using Test

const MANUAL_BISECTION_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const MANUAL_BISECTION_JOB = joinpath(
    MANUAL_BISECTION_ROOT, "scripts", "analysis", "pnjl_c2_cep_manual_bisection_job.jl",
)

@testset "C2 manual bisection cost validation contract" begin
    source = read(MANUAL_BISECTION_JOB, String)
    @test occursin("_validate_manual_cost_rows(cost_rows)", source)
    @test occursin("all(item -> item.failed_points == 0, cost_rows)", source)
    @test occursin(
        "all(item -> item.point_requests == item.unique_solves + item.cache_hits, cost_rows)",
        source,
    )
    @test !occursin("all(item -> item.failed_points == 0 for item in cost_rows)", source)
    @test !occursin(
        "all(item -> item.point_requests == item.unique_solves + item.cache_hits for item in cost_rows)",
        source,
    )
end
