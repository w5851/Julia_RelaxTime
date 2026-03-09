using Test

if !isdefined(Main, :_compute_relaxtime_tau_point)
    include(joinpath(@__DIR__, "common", "tau_validation_helpers.jl"))
end

const RELAXTIME_VALIDATION_DATA_PATH = validation_targets_path(
    "relaxtime",
    "literature",
    "tau",
    "relaxtime_tau_literature_targets_v1.csv",
)

@testset "RelaxTime literature digitized tau targets" begin
    targets = _load_tau_targets(RELAXTIME_VALIDATION_DATA_PATH)
    @test !isempty(targets)

    for row in targets
        tau = _compute_relaxtime_tau_point(row.T_MeV, row.muB_MeV)
        species = _tau_species_key(row)
        actual_tau_fm = getproperty(tau, species)

        @test isfinite(actual_tau_fm)
        @test isapprox(actual_tau_fm, row.expected_tau_fm; rtol=row.rtol, atol=0.0)
    end
end