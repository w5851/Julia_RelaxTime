using Test

if !isdefined(Main, :_load_legacy_transport_scalar_targets)
    include(joinpath(@__DIR__, "legacy_transport_validation_helpers.jl"))
end

const RELAXTIME_ZETA_OVER_S_LEGACY_GUARDRAIL_PATH = validation_targets_path(
    "relaxtime",
    "legacy",
    "transport",
    "relaxtime_zeta_over_s_legacy_targets_v1.csv",
)

@testset "RelaxTime legacy zeta_over_s guardrail" begin
    targets = _load_legacy_transport_scalar_targets(RELAXTIME_ZETA_OVER_S_LEGACY_GUARDRAIL_PATH)
    @test !isempty(targets)

    for row in targets
        @test row.field == "zeta_over_s"
        transport_point = _compute_relaxtime_literature_transport_point_bulk(row.T_MeV, row.muB_MeV)
        actual = transport_point.transport.zeta / transport_point.s_fm3inv
        @test isfinite(actual)
        @test isapprox(actual, row.expected_value; rtol=row.rtol, atol=0.0)
    end
end