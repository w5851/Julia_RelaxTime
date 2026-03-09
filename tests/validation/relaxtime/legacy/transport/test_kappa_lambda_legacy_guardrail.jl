using Test

if !isdefined(Main, :_load_legacy_transport_scalar_targets)
    include(joinpath(@__DIR__, "legacy_transport_validation_helpers.jl"))
end

const RELAXTIME_KAPPA_LAMBDA_LEGACY_GUARDRAIL_PATH = validation_targets_path(
    "relaxtime",
    "legacy",
    "transport",
    "relaxtime_kappa_lambda_fixedpoint_targets_v1.csv",
)

@testset "RelaxTime legacy kappa/lambda guardrail" begin
    targets = _load_legacy_transport_scalar_targets(RELAXTIME_KAPPA_LAMBDA_LEGACY_GUARDRAIL_PATH)
    @test !isempty(targets)

    for row in targets
        @test row.field in ("kappa_BB", "kappa_QQ", "kappa_SS", "lambda")
        actual_transport = _compute_relaxtime_literature_transport_point(row.T_MeV, row.muB_MeV).transport
        actual = getproperty(actual_transport, Symbol(row.field))
        @test isfinite(actual)
        @test isapprox(actual, row.expected_value; rtol=row.rtol, atol=0.0)
    end
end