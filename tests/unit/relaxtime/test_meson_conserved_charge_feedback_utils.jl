using Test

const _MESON_FEEDBACK_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(
    _MESON_FEEDBACK_ROOT,
    "scripts",
    "analysis",
    "relaxtime",
    "meson_conserved_charge_feedback_utils.jl",
))
using .MesonConservedChargeFeedbackUtils

@testset "charged meson conserved-charge feedback helpers" begin
    chemical = charged_meson_chemical_potentials(-0.05, 0.4)
    @test chemical.mu_pi_plus == -0.05
    @test chemical.mu_pi_minus == 0.05
    @test chemical.mu_K_plus ≈ 0.35
    @test chemical.mu_K_minus ≈ -0.35
    @test_throws ArgumentError charged_meson_chemical_potentials(NaN, 0.0)

    meson = charged_meson_conserved_densities(0.3, 0.2, 0.1, 0.04)
    @test meson.rho_B == 0.0
    @test meson.rho_Q ≈ 0.16
    @test meson.rho_S ≈ 0.06
    @test_throws ArgumentError charged_meson_conserved_densities(-1e-3, 0.0, 0.0, 0.0)

    residual = total_conserved_charge_residual(
        (rho_B=0.1, rho_Q=0.02, rho_S=-0.01),
        (rho_Q=0.02, rho_S=0.01);
        charge_to_baryon_ratio=0.4,
        rho0=0.16,
    )
    @test residual.charge_raw ≈ 0.0 atol=1e-15
    @test residual.strangeness_raw ≈ 0.0 atol=1e-15
    @test residual.norm ≈ 0.0 atol=1e-15
    @test_throws ArgumentError total_conserved_charge_residual(
        (rho_B=0.1, rho_Q=0.02, rho_S=0.0),
        (rho_Q=0.0, rho_S=0.0);
        rho0=0.0,
    )

    target_mu_Q = -0.08
    target_mu_S = 0.31
    linear_evaluator = (mu_Q, mu_S) -> (
        rho_B_q=0.12,
        rho_Q_q=0.4 * 0.12 + 0.16 * (mu_Q - target_mu_Q),
        rho_S_q=0.16 * (mu_S - target_mu_S),
        n_pi_plus=0.0,
        n_pi_minus=0.0,
        n_K_plus=0.0,
        n_K_minus=0.0,
    )
    solved = solve_outer_conserved_charge_feedback(
        linear_evaluator,
        0.0,
        0.0;
        residual_tolerance=1e-10,
        finite_difference_step=1e-4,
        maximum_step=1.0,
    )
    @test solved.converged
    @test solved.reason === :residual_tolerance
    @test solved.mu_Q ≈ target_mu_Q atol=1e-9
    @test solved.mu_S ≈ target_mu_S atol=1e-9
    @test solved.residual_norm <= 1e-10
    @test solved.evaluation_count <= 4

    singular_evaluator = (mu_Q, mu_S) -> (
        rho_B_q=0.1,
        rho_Q_q=0.0,
        rho_S_q=0.1,
        n_pi_plus=0.0,
        n_pi_minus=0.0,
        n_K_plus=0.0,
        n_K_minus=0.0,
    )
    singular = solve_outer_conserved_charge_feedback(singular_evaluator, 0.0, 0.0)
    @test !singular.converged
    @test singular.reason === :singular_jacobian

    invalid_evaluator = (mu_Q, mu_S) -> (
        rho_B_q=0.1,
        rho_Q_q=0.04,
        rho_S_q=0.0,
        n_pi_plus=-0.1,
        n_pi_minus=0.0,
        n_K_plus=0.0,
        n_K_minus=0.0,
    )
    invalid = solve_outer_conserved_charge_feedback(invalid_evaluator, 0.0, 0.0)
    @test !invalid.converged
    @test invalid.reason === :initial_evaluation_failed
    @test occursin("n_pi_plus", invalid.message)

    @test_throws ArgumentError solve_outer_conserved_charge_feedback(
        linear_evaluator,
        0.0,
        0.0;
        max_evaluations=2,
    )
end
