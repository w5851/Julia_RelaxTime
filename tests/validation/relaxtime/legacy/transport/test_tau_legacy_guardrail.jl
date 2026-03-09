using Test

if !isdefined(Main, :_compute_relaxtime_tau_point)
    include(joinpath(@__DIR__, "legacy_transport_validation_helpers.jl"))
end

const RELAXTIME_TAU_LEGACY_GUARDRAIL_PATH = validation_targets_path(
    "relaxtime",
    "legacy",
    "transport",
    "relaxtime_tau_legacy_targets_v1.csv",
)

@testset "RelaxTime legacy tau guardrail" begin
    targets = _load_tau_targets(RELAXTIME_TAU_LEGACY_GUARDRAIL_PATH)
    @test !isempty(targets)

    for row in targets
        tau = _compute_relaxtime_tau_point(row.T_MeV, row.muB_MeV)
        species = _tau_species_key(row)
        actual_tau_fm = getproperty(tau, species)

        @test isfinite(actual_tau_fm)
        @test isapprox(actual_tau_fm, row.expected_tau_fm; rtol=row.rtol, atol=0.0)
    end
end