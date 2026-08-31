using Test

const _CHARGED_RPA_NEGATIVE_DENSITY_SCRIPT = normpath(joinpath(
    @__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime",
    "diagnose_charged_rpa_bu_negative_density.jl",
))

@testset "charged-RPA/BU negative-density diagnostic contract" begin
    @test isfile(_CHARGED_RPA_NEGATIVE_DENSITY_SCRIPT)
    source = read(_CHARGED_RPA_NEGATIVE_DENSITY_SCRIPT, String)
    @test Meta.parseall(source) isa Expr
    @test occursin("FixedMuBConservedCharges", source)
    @test occursin("stable_meson_number_density", source)
    @test occursin("strict_bw_meson_number_density", source)
    @test occursin("strict_bw_qpole_meson_number_density", source)
    @test occursin("phase_shift_meson_number_density", source)
    @test occursin("bose_support_gate", source)
    @test occursin("phase_shift_gbu_reference", source)
    @test occursin("arg_propagator", source)
    @test occursin("arg_inverse_propagator", source)
    @test occursin("high_energy_zero", source)
    @test occursin("legacy_domega_over_2pi", source)
    @test occursin("phase_shell_breakdown", source)
    @test occursin("accepted_flags", source)
    @test occursin("CHARGED_RPA_NEGATIVE_DENSITY_OUTPUT", source)
    @test !occursin("production_candidate_status=\"authorized\"", source)
end
