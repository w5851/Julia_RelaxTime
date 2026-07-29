using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Rho-support cascade router" begin
    @test Models.RhoSupportConfig().target_point_count == 9
    @test Models.RhoSupportConfig().max_extra_points == 12
    @test_throws ArgumentError Models.analyze_rho_support_cascade(
        collect(0.0:0.25:3.0), collect(0.0:0.25:3.0);
        config=Models.RhoSupportConfig(target_point_count=4),
    )

    rho = collect(0.0:0.25:3.0)
    monotone = Models.analyze_rho_support_cascade(rho, rho)
    @test monotone.status in (:unresolved, :near_critical)
    @test monotone.extra_point_count == 0

    x = rho .- 1.5
    s_curve = x .^ 3 .- 0.45 .* x
    routed = Models.analyze_rho_support_cascade(
        rho,
        s_curve;
        sample_mu = value -> value^3 - 0.45 * (value - 1.5),
    )
    @test routed.total_point_count >= routed.coarse_point_count
    @test routed.extra_point_count <= Models.RhoSupportConfig().max_extra_points

    prior = Models.RhoSupportPrior(100.0, 1.5, 0.75)
    prior_routed = Models.analyze_rho_support_cascade(
        rho, monotone.baseline_has_s_shape ? s_curve : rho;
        prior=prior,
        sample_mu=value -> value,
    )
    @test prior_routed.support_origin in (:temperature_prior, :coarse_and_prior, :coarse_low_slope, :none)
end
