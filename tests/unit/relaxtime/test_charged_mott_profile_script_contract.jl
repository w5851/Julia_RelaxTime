using Test

const _CHARGED_MOTT_SCRIPT = joinpath(
    @__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime",
    "audit_charged_mott_profiles.jl",
)

@testset "charged Mott profile audit contract" begin
    @test isfile(_CHARGED_MOTT_SCRIPT)
    source = read(_CHARGED_MOTT_SCRIPT, String)
    @test Meta.parseall(source) isa Expr
    @test occursin("strict_mott_gate", source)
    @test occursin("before_bound_state_count", source)
    @test occursin("after_bound_state_count", source)
    @test occursin("before_mott_gap_inv_fm", source)
    @test occursin("bound_state_q_count_min", source)
    @test occursin("bound_state_complex_q_count", source)
    @test occursin("mott_gate_passed", source)
    @test occursin("production_candidate_status=\"not_authorized\"", source)
end
