using Test

if !isdefined(Main, :_load_legacy_transport_scalar_targets)
    include(joinpath(@__DIR__, "legacy_transport_validation_helpers.jl"))
end

const RELAXTIME_SIGMA_T_LEGACY_GUARDRAIL_PATH = validation_targets_path(
    "relaxtime",
    "legacy",
    "transport",
    "relaxtime_sigma_t_legacy_targets_v1.csv",
)

@testset "RelaxTime legacy sigma_t guardrail" begin
    targets = _load_legacy_transport_scalar_targets(RELAXTIME_SIGMA_T_LEGACY_GUARDRAIL_PATH)
    @test !isempty(targets)

    for row in targets
        @test row.field == "sigma_t"
        actual = _compute_relaxtime_literature_transport_point(row.T_MeV, row.muB_MeV).transport.sigma / row.T_MeV * Main.Constants_PNJL.ħc_MeV_fm
        @test isfinite(actual)
        @test isapprox(actual, row.expected_value; rtol=row.rtol, atol=0.0)
    end
end