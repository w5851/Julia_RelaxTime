using Test

const _CHARGED_RPA_FIXEDPOINT_SCRIPT = normpath(joinpath(
    @__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime",
    "compare_charged_rpa_ordered_fixedpoints.jl",
))

@testset "charged-RPA ordered fixed-point script contract" begin
    @test isfile(_CHARGED_RPA_FIXEDPOINT_SCRIPT)
    source = read(_CHARGED_RPA_FIXEDPOINT_SCRIPT, String)
    @test Meta.parseall(source) !== nothing
    @test occursin("FixedMuBConservedCharges", source)
    @test occursin("production_candidate_status=\"not_authorized\"", source)
    @test occursin("charge_symmetric_reference", source)
    @test occursin("finite_bqs_reference", source)
    @test occursin("zero_density_degenerate_reference", source)
    @test occursin("(:pi_plus, :pi_minus, :K_plus, :K_minus)", source)
    @test occursin("ordered_retarded_eta_coarse", source)
    @test occursin("ordered_retarded_nodes_coarse", source)
    @test occursin("ordered_retarded_refined", source)
    @test occursin("ordered_legacy_B0", source)
    @test occursin("legacy_symmetrized_B0", source)
    @test occursin("charged_rpa_coupling(kernel, spec)", source)
    @test occursin("kernel_pair=String(spec.kernel_pair)", source)
    @test occursin("Pi_relative_to_refined", source)
    @test occursin("inverse_relative_to_refined", source)
    @test occursin("propagator_relative_to_refined", source)
    @test occursin("Pi_positive_energy_partner_relative", source)
    @test occursin("propagator_positive_energy_partner_relative", source)
    @test occursin("not_applicable_finite_mu", source)
end
