using Test

const _PHASE_GUIDED_SCAN_SCRIPT = joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_phase_guided_transport_scan.jl")
const _PHASE_GUIDED_PRODUCTION_WORKFLOW = joinpath(@__DIR__, "..", "..", "..", ".github", "workflows", "relaxtime-phase-guided-transport-production.yml")
if !isdefined(Main, :run_phase_guided_scan)
    include(_PHASE_GUIDED_SCAN_SCRIPT)
end

@testset "phase-guided transport plan mode a" begin
    opts = Main.PhaseGuidedTransportScanCLI.PhaseGuidedScanOptions(
        :mode_a_fixed_muB_phase_scaled,
        "tmp",
        "unit_case_a",
        [-0.5, -0.2, 0.0, 0.2],
        [0.0, 450.0, 900.0],
        [1.0, 1.1],
        Float64[],
        :match_thermo,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        false,
        true,
        true,
        false,
        false,
    )
    plan = Main.PhaseGuidedTransportScanPlan.build_plan(opts)
    @test plan.total == 24
    @test unique(p.plot_panel for p in plan.points if p.muB_MeV == 0.0) == ["muB0.0"]
    @test unique(p.plot_series for p in plan.points if p.alpha_T == 1.0) == ["alpha1.0"]
    @test any(occursin("T=", p.plot_series_label) for p in plan.points if p.muB_MeV == 450.0)
    refs_0 = unique(p.phase_reference_kind for p in plan.points if p.muB_MeV == 0.0)
    refs_450 = unique(p.phase_reference_kind for p in plan.points if p.muB_MeV == 450.0)
    refs_900 = unique(p.phase_reference_kind for p in plan.points if p.muB_MeV == 900.0)
    @test refs_0 == [:crossover]
    @test refs_450 == [:crossover]
    @test refs_900 == [:first_order]
    @test all(isfinite(p.T_phase_base_MeV) for p in plan.points)
    base_450 = unique(p.T_phase_base_MeV for p in plan.points if p.muB_MeV == 450.0)
    @test length(base_450) == 1
end

@testset "phase-guided transport plan mode b" begin
    opts = Main.PhaseGuidedTransportScanCLI.PhaseGuidedScanOptions(
        :mode_b_fixed_T_sparse_muB,
        "tmp",
        "unit_case_b",
        [-0.5, -0.2, 0.0, 0.2],
        [0.0, 900.0],
        Float64[1.0],
        [120.0, 130.0, 138.0],
        :match_thermo,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        false,
        true,
        true,
        false,
        false,
    )
    plan = Main.PhaseGuidedTransportScanPlan.build_plan(opts)
    @test plan.total == 24
    @test unique(p.plot_panel for p in plan.points if p.T_MeV == 120.0) == ["T120.0"]
    @test unique(p.plot_series for p in plan.points if p.muB_MeV == 900.0) == ["muB900.0"]
    @test any(p.phase_reference_kind == :cep_neighbor for p in plan.points)
    @test any(p.phase_reference_kind == :first_order for p in plan.points)
    @test any(p.phase_reference_kind == :crossover for p in plan.points)
end

@testset "phase-guided transport crossover reference uses muq" begin
    xi = 0.0
    muq_ref_MeV = 90.0
    muB_MeV = 3.0 * muq_ref_MeV

    direct_T, direct_xi = Main.GapTransportScanPhaseEquilibrium.interpolate_crossover_temperature(xi, muq_ref_MeV)
    via_plan_T, via_plan_xi = Main.PhaseGuidedTransportScanPlan._interpolate_crossover_temperature(muB_MeV, xi)

    @test isfinite(direct_T)
    @test isfinite(via_plan_T)
    @test direct_xi == via_plan_xi
    @test isapprox(via_plan_T, direct_T; atol=1e-12, rtol=0.0)
end

@testset "phase-guided production workflow forwards required CLI args" begin
    workflow_text = read(_PHASE_GUIDED_PRODUCTION_WORKFLOW, String)

    @test occursin("--case-name", workflow_text)
    @test occursin("canonical_xi_list=\"-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5\"", workflow_text)
    @test occursin("--xi-list \"\$canonical_xi_list\"", workflow_text)
    @test occursin("--propagator-xi-policy", workflow_text)
    @test occursin("--channel-diagnostics", workflow_text)
    @test occursin("--sigma-grid-n", workflow_text)
end
