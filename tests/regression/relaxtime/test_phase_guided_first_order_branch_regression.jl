using Test

const PROJECT_ROOT_PHASE_GUIDED_BRANCH = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCAN_SCRIPT_PHASE_GUIDED_BRANCH = joinpath(
    PROJECT_ROOT_PHASE_GUIDED_BRANCH,
    "scripts",
    "relaxtime",
    "run_gap_transport_scan.jl",
)

if !isdefined(Main, :execute_gap_transport_scan_point!)
    include(SCAN_SCRIPT_PHASE_GUIDED_BRANCH)
end

@testset "phase-guided first-order points select the stable branch" begin
    T_MeV = 125.06991551
    muB_MeV = 900.0
    xi = -0.01
    T_fm = T_MeV / Main.ħc_MeV_fm
    muq_fm = (muB_MeV / 3.0) / Main.ħc_MeV_fm

    opts = Main.parse_args([
        "--output", "unused.csv",
        "--tmin", string(T_MeV),
        "--tmax", string(T_MeV),
        "--mubmin", string(muB_MeV),
        "--mubmax", string(muB_MeV),
        "--xi-list", string(xi),
        "--p-num", "12",
        "--t-num", "6",
    ])

    continuation = Main.GapTransportScanPhaseEquilibrium.solve_models_equilibrium(
        T_fm,
        muq_fm,
        xi,
        Main.Models.QUARK_SEED_5,
        opts,
    )
    @test continuation !== nothing
    @test Main.GapTransportScanPhaseEquilibrium.classify_phase_from_solution(continuation) === :quark

    selected, _, selected_phase, diag = Main.GapTransportScanPhaseEquilibrium.solve_equilibrium_with_diagnostics(
        T_MeV,
        muB_MeV,
        xi,
        opts;
        previous_solution=continuation.x_state,
        previous_phase=:quark,
    )

    @test selected_phase === :hadron
    @test Main.GapTransportScanPhaseEquilibrium.classify_phase_from_solution(selected) === :hadron
    @test selected.omega < continuation.omega
    @test selected.masses[1] > continuation.masses[1]
    @test diag.phase_structure === :first_order_possible
    @test diag.phase_curr_hint === :hadron
    @test diag.equilibrium_selection_policy === :pressure_max_under_constraints
    @test startswith(diag.seed_source, "first_order_stable_multiseed_")

    bulk = Main.Models.bulk_viscosity_coefficients(
        T_fm,
        muq_fm;
        xi=xi,
        p_num=12,
        t_num=6,
        model=Main.PNJL_MODEL,
        base_state=selected.x_state,
        base_masses=selected.masses,
        base_mu_vec=selected.mu_vec,
        base_state_source=:workflow_equilibrium,
    )
    @test bulk.base_state_source === :workflow_equilibrium
    @test bulk.primal_solve_count == 0
    @test bulk.branch_locked
    @test all(isapprox.(bulk.masses, selected.masses; rtol=1e-4, atol=1e-8))

    anchor = Main.GapTransportScanPhaseEquilibrium.direct_coexistence_anchor(
        muB_MeV,
        T_MeV,
        opts,
    )
    @test isapprox(anchor.T_mid_MeV, 125.76660544; atol=2e-4, rtol=0.0)
    @test anchor.bracket_width_MeV < 1e-5
    @test anchor.delta_omega_lower <= 0.0 <= anchor.delta_omega_upper

    certification = Main.GapTransportScanPhaseEquilibrium.certify_coexistence_side_points(
        anchor,
        muB_MeV,
        opts,
    )
    @test certification.certified
    @test certification.delta_xi == 0.003
    @test certification.anchor_convergence_certified
    @test abs(certification.anchor_convergence_delta_MeV) < 0.1
    @test certification.node_configs == [(p_num=12, t_num=6), (p_num=24, t_num=8)]
    @test certification.convergence_anchor.p_num == 24
    @test certification.convergence_anchor.t_num == 8
    @test all(row.minus_certified && row.plus_certified for row in certification.evidence)
end
