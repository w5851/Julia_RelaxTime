using Test

const _MESON_FEEDBACK_ANALYSIS_ROOT = normpath(joinpath(
    @__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime",
))

@testset "meson conserved-charge outer-feedback analysis contract" begin
    helper_path = joinpath(_MESON_FEEDBACK_ANALYSIS_ROOT, "meson_conserved_charge_feedback_utils.jl")
    script_path = joinpath(_MESON_FEEDBACK_ANALYSIS_ROOT, "meson_conserved_charge_outer_feedback_spike.jl")
    @test isfile(helper_path)
    @test isfile(script_path)

    helper_source = read(helper_path, String)
    script_source = read(script_path, String)
    @test Meta.parseall(helper_source) !== nothing
    @test Meta.parseall(script_source) !== nothing
    @test occursin("FixedMuBConservedCharges", script_source)
    @test occursin("flavor_mu_from_bqs", script_source)
    @test occursin("model_rho", script_source)
    @test occursin("scheme=:current", script_source)
    @test occursin("density_policy=:x_min_cut", script_source)
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
