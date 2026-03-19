using Test

const _SCAN_QUALITY_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "scan_quality.jl"))
Base.include(Main, _SCAN_QUALITY_PATH)

using Main.ScanQuality: assess_point_quality

@testset "scan quality assessment" begin
    @testset "healthy point is not flagged" begin
        tau = (u=1.0, d=1.1, s=1.2, ubar=1.05, dbar=1.15, sbar=1.25)
        flagged, reason, metric = assess_point_quality(tau, 0.25, 0.4)
        @test flagged == false
        @test reason == "ok"
        @test isfinite(metric)
    end

    @testset "nonpositive tau is flagged" begin
        tau = (u=0.0, d=1.1, s=1.2, ubar=1.05, dbar=1.15, sbar=1.25)
        flagged, reason, metric = assess_point_quality(tau, 0.25, 0.4)
        @test flagged == true
        @test reason == "tau_invalid_u"
        @test metric == 0.0
    end

    @testset "large tau_u/tau_ubar ratio is flagged" begin
        tau = (u=12.0, d=1.1, s=1.2, ubar=1.5, dbar=1.15, sbar=1.25)
        flagged, reason, metric = assess_point_quality(tau, 0.25, 0.4)
        @test flagged == true
        @test reason == "tau_u_ubar_ratio_high"
        @test metric > 6.0
    end

    @testset "nonfinite eta_over_s is flagged" begin
        tau = (u=1.0, d=1.1, s=1.2, ubar=1.05, dbar=1.15, sbar=1.25)
        flagged, reason, metric = assess_point_quality(tau, NaN, 0.4)
        @test flagged == true
        @test reason == "eta_over_s_nonfinite"
        @test isnan(metric)
    end
end
