using Test

const _REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _SCRIPT_PATH = joinpath(_REPO_ROOT, "scripts", "relaxtime", "evaluate_xi_smoothness.jl")

@testset "xi smoothness metrics" begin
    @test isfile(_SCRIPT_PATH)

    if isfile(_SCRIPT_PATH)
        Base.include(Main, _SCRIPT_PATH)

        @testset "S2 is zero on constant sequence" begin
            ys = [1.0, 1.0, 1.0, 1.0, 1.0]
            @test isapprox(compute_s2(ys), 0.0; atol=1e-12, rtol=0.0)
        end

        @testset "S2 is positive on spike sequence" begin
            ys = [1.0, 1.0, 5.0, 1.0, 1.0]
            @test compute_s2(ys) > 0.0
        end

        @testset "N_spike counting" begin
            ys = [1.0, 1.0, 5.0, 1.0, 1.0]
            @test count_spikes(ys; spike_threshold=1.0) == 3
            @test count_spikes(fill(2.0, 5); spike_threshold=1.0) == 0
        end

        @testset "classify dual-threshold logic" begin
            thresholds = (
                smooth_s2 = 0.02,
                suspect_s2 = 0.10,
                smooth_s1jump = 0.10,
                suspect_s1jump = 0.30,
                suspect_max_spikes = 2,
            )

            smooth_label, _ = classify_smoothness(0.01, 0.05, 0;
                smooth_s2=thresholds.smooth_s2,
                suspect_s2=thresholds.suspect_s2,
                smooth_s1jump=thresholds.smooth_s1jump,
                suspect_s1jump=thresholds.suspect_s1jump,
                suspect_max_spikes=thresholds.suspect_max_spikes,
            )
            suspect_label, _ = classify_smoothness(0.06, 0.20, 1;
                smooth_s2=thresholds.smooth_s2,
                suspect_s2=thresholds.suspect_s2,
                smooth_s1jump=thresholds.smooth_s1jump,
                suspect_s1jump=thresholds.suspect_s1jump,
                suspect_max_spikes=thresholds.suspect_max_spikes,
            )
            not_smooth_label, _ = classify_smoothness(0.12, 0.20, 1;
                smooth_s2=thresholds.smooth_s2,
                suspect_s2=thresholds.suspect_s2,
                smooth_s1jump=thresholds.smooth_s1jump,
                suspect_s1jump=thresholds.suspect_s1jump,
                suspect_max_spikes=thresholds.suspect_max_spikes,
            )

            @test smooth_label == "smooth"
            @test suspect_label == "suspect"
            @test not_smooth_label == "not_smooth"
        end
    end
end
