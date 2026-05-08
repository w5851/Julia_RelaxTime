using Test

const _PHASE_GUIDED_SCAN_SCRIPT = joinpath(@__DIR__, "..", "..", "..", "scripts", "relaxtime", "run_phase_guided_transport_scan.jl")
if !isdefined(Main, :run_phase_guided_scan)
    include(_PHASE_GUIDED_SCAN_SCRIPT)
end

@testset "phase-guided transport plan mode a" begin
    opts = Main.PhaseGuidedTransportScanCLI.PhaseGuidedScanOptions(
        :mode_a_fixed_muB_phase_scaled,
        "tmp",
        "unit_case_a",
        [-0.5, -0.2, 0.0, 0.2],
        [0.0, 400.0],
        [1.0, 1.1],
        Float64[],
        true,
        true,
        false,
        false,
    )
    plan = Main.PhaseGuidedTransportScanPlan.build_plan(opts)
    @test plan.total == 16
    refs_0 = unique(p.phase_reference_kind for p in plan.points if p.muB_MeV == 0.0)
    refs_400 = unique(p.phase_reference_kind for p in plan.points if p.muB_MeV == 400.0)
    # With crossover_dense.csv present for the current xi anchors, mu_B = 0
    # now consumes crossover references directly.
    @test refs_0 == [:crossover]
    @test :first_order in refs_400
    @test all(isfinite(p.T_phase_base_MeV) for p in plan.points)
end

@testset "phase-guided transport plan mode b" begin
    opts = Main.PhaseGuidedTransportScanCLI.PhaseGuidedScanOptions(
        :mode_b_fixed_T_sparse_muB,
        "tmp",
        "unit_case_b",
        [-0.5, -0.2, 0.0, 0.2],
        [0.0, 260.0, 420.0],
        Float64[1.0],
        [120.0, 138.0],
        true,
        true,
        false,
        false,
    )
    plan = Main.PhaseGuidedTransportScanPlan.build_plan(opts)
    @test plan.total == 24
    @test any(p.phase_reference_kind == :cep_neighbor for p in plan.points)
    @test any(p.phase_reference_kind == :first_order for p in plan.points)
    @test any(p.phase_reference_kind == :crossover for p in plan.points)
end
