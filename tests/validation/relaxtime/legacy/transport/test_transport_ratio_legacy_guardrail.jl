using Test

if !isdefined(Main, :_load_legacy_transport_scalar_targets)
    include(joinpath(@__DIR__, "legacy_transport_validation_helpers.jl"))
end

const RELAXTIME_TRANSPORT_RATIO_LEGACY_GUARDRAIL_PATH = validation_targets_path(
    "relaxtime",
    "legacy",
    "transport",
    "relaxtime_transport_ratio_legacy_targets_v1.csv",
)

@testset "RelaxTime legacy transport ratio guardrail" begin
    targets = _load_legacy_transport_scalar_targets(RELAXTIME_TRANSPORT_RATIO_LEGACY_GUARDRAIL_PATH)
    @test !isempty(targets)

    for row in targets
        actual = if row.field == "viscous_conductive_coupling_ratio"
            _compute_relaxtime_literature_transport_point(row.T_MeV, row.muB_MeV).transport.viscous_conductive_coupling_ratio
        elseif row.field == "bulk_to_shear_viscosity_ratio"
            _compute_relaxtime_literature_transport_point_bulk(row.T_MeV, row.muB_MeV).transport.bulk_to_shear_viscosity_ratio
        else
            error("unexpected transport ratio field: $(row.field)")
        end

        @test isfinite(actual)
        @test isapprox(actual, row.expected_value; rtol=row.rtol, atol=0.0)
    end
end