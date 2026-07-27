using Test

const _MESON_FEEDBACK_ANALYSIS_ROOT = normpath(joinpath(
    @__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime",
))

@testset "meson conserved-charge outer-feedback analysis contract" begin
    helper_path = joinpath(_MESON_FEEDBACK_ANALYSIS_ROOT, "meson_conserved_charge_feedback_utils.jl")
    runtime_path = joinpath(_MESON_FEEDBACK_ANALYSIS_ROOT, "meson_conserved_charge_feedback_runtime.jl")
    script_path = joinpath(_MESON_FEEDBACK_ANALYSIS_ROOT, "meson_conserved_charge_outer_feedback_spike.jl")
    scan_path = joinpath(_MESON_FEEDBACK_ANALYSIS_ROOT, "scan_meson_conserved_charge_feedback_freezeout.jl")
    benchmark_path = normpath(joinpath(_MESON_FEEDBACK_ANALYSIS_ROOT, "..", "..", "perf", "relaxtime", "bench_meson_conserved_charge_partial_feedback.jl"))
    @test isfile(helper_path)
    @test isfile(runtime_path)
    @test isfile(script_path)
    @test isfile(scan_path)
    @test isfile(benchmark_path)

    helper_source = read(helper_path, String)
    runtime_source = read(runtime_path, String)
    script_source = read(script_path, String)
    scan_source = read(scan_path, String)
    benchmark_source = read(benchmark_path, String)
    @test Meta.parseall(helper_source) !== nothing
    @test Meta.parseall(runtime_source) !== nothing
    @test Meta.parseall(script_source) !== nothing
    @test Meta.parseall(scan_source) !== nothing
    @test Meta.parseall(benchmark_source) !== nothing
    @test occursin("FixedMuBConservedCharges", runtime_source)
    @test occursin("flavor_mu_from_bqs", runtime_source)
    @test occursin("model_rho", runtime_source)
    @test occursin("scheme=:current", runtime_source)
    @test occursin("density_policy=settings.density_policy", runtime_source)
    @test occursin("FeedbackSettings", runtime_source)
    @test occursin("candidate_timing_summary", runtime_source)
    @test occursin("meson_conserved_charge_feedback_runtime.jl", script_source)
    @test occursin("choose_freezeout_sqrts_grid", scan_source)
    @test occursin("baseline_freezeout", scan_source)
    @test occursin("Kplus_over_piplus", scan_source)
    @test occursin("Kminus_over_piminus", scan_source)
    @test occursin("MESON_FEEDBACK_BENCH_SUMMARY", scan_source)
    @test occursin("fresh evaluator/cache", benchmark_source)
    @test occursin("warmup", lowercase(benchmark_source))
    @test occursin("candidate_elapsed_s", benchmark_source)
    @test occursin("recommended_sqrts_grid_GeV", benchmark_source)
    @test occursin("MESON_FEEDBACK_RUN_REFINED_OUTER", script_source)
    @test occursin("MESON_FEEDBACK_INITIAL_MUQ_MEV", script_source)

    for field in (
        "rho_Q_M",
        "rho_S_M",
        "rho_Q_total",
        "rho_S_total",
        "n_pi_plus",
        "n_pi_minus",
        "n_K_plus",
        "n_K_minus",
        "Kplus_over_piplus",
        "Kminus_over_piminus",
        "pi_plus_min_E_minus_mu",
        "K_plus_min_E_minus_mu",
        "gap_residual_norm",
        "outer_initial_muQ_MeV",
        "unique_candidate_count",
        "verdict",
    )
        @test occursin(field, script_source)
    end
end
