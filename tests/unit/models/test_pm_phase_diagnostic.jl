using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "PM phase diagnostic seed/type contracts" begin
    @test isdefined(Models, :PM_BRANCH_STATUSES)
    @test Models.PM_BRANCH_STATUSES == (:accepted, :nonconverged, :branch_jump, :invalid_thermo)

    @test isdefined(Models, :PM_SEED_SOURCES)
    @test Models.PM_SEED_SOURCES == (:seed0, :previous_same_branch, :manual_override)

    @test isdefined(Models, :PM_ENDPOINT_CAUSES)
    @test Models.PM_ENDPOINT_CAUSES == (:physical_loss_candidate, :nonconvergence, :branch_jump, :out_of_grid)

    @test isdefined(Models, :PM_COMPARISON_STATUSES)
    @test Models.PM_COMPARISON_STATUSES == (:both, :pm_only, :maxwell_only, :neither)

    @test isdefined(Models, :normalize_pm_seed_pair)
    seed_pair = Models.normalize_pm_seed_pair((
        hadron_seed0=[1.0, 2.0, 3.0],
        quark_seed0=[4.0, 5.0, 6.0],
    ))
    @test seed_pair.hadron_seed0 == [1.0, 2.0, 3.0]
    @test seed_pair.quark_seed0 == [4.0, 5.0, 6.0]
    @test seed_pair.continuity_mode == :branch_local
    @test seed_pair.fallback_mode == :none

    struct_seed_pair = Models.PMSeedPair(
        hadron_seed0=[7.0, 8.0],
        quark_seed0=[9.0, 10.0],
        continuity_mode=:branch_local,
        fallback_mode=:none,
    )
    normalized_struct = Models.normalize_pm_seed_pair(struct_seed_pair)
    @test normalized_struct.hadron_seed0 == [7.0, 8.0]
    @test normalized_struct.quark_seed0 == [9.0, 10.0]

    @test isdefined(Models, :pm_next_seed_source)
    @test Models.pm_next_seed_source(true, :hadron) == :previous_same_branch
    @test Models.pm_next_seed_source(false, :hadron) == :seed0
end

@testset "PM phase diagnostic summary helpers" begin
    @test isdefined(Models, :_pm_interpolate_transition_mu)
    rows_transition = [
        (mu_MeV=291.10, delta_omega=0.12, delta_pressure=-0.12, hadron_exists=true, quark_exists=true, hadron_status=:accepted, quark_status=:accepted),
        (mu_MeV=291.20, delta_omega=-0.08, delta_pressure=0.08, hadron_exists=true, quark_exists=true, hadron_status=:accepted, quark_status=:accepted),
    ]
    mu_transition = Models._pm_interpolate_transition_mu(rows_transition)
    @test isapprox(mu_transition, 291.16; atol=1e-8)

    @test isdefined(Models, :_pm_extract_endpoints)
    rows_window = [
        (mu_MeV=290.9, hadron_exists=true, quark_exists=false, hadron_status=:accepted, quark_status=:out_of_grid),
        (mu_MeV=291.0, hadron_exists=true, quark_exists=true, hadron_status=:accepted, quark_status=:accepted),
        (mu_MeV=291.1, hadron_exists=true, quark_exists=true, hadron_status=:accepted, quark_status=:accepted),
        (mu_MeV=291.2, hadron_exists=false, quark_exists=true, hadron_status=:branch_jump, quark_status=:accepted),
    ]
    endpoints = Models._pm_extract_endpoints(rows_window)
    @test endpoints.hadron_endpoint_mu_MeV == 291.1
    @test endpoints.quark_endpoint_mu_MeV == 291.0
    @test isapprox(endpoints.bistable_window_width_MeV, 0.1; atol=1e-10)
    @test endpoints.branch_disappears_first == :hadron

    @test isdefined(Models, :_pm_has_bistable_window)
    @test Models._pm_has_bistable_window(rows_window)
    rows_nowindow = [
        (mu_MeV=290.9, hadron_exists=true, quark_exists=false, hadron_status=:accepted, quark_status=:out_of_grid),
        (mu_MeV=291.0, hadron_exists=true, quark_exists=false, hadron_status=:accepted, quark_status=:out_of_grid),
        (mu_MeV=291.1, hadron_exists=false, quark_exists=true, hadron_status=:nonconverged, quark_status=:accepted),
    ]
    @test !Models._pm_has_bistable_window(rows_nowindow)

    @test isdefined(Models, :_pm_compare_with_maxwell)
    cmp = Models._pm_compare_with_maxwell(291.16, 291.19; comparison_mu_tol=0.05)
    @test cmp.comparison_status == :both
    @test isapprox(cmp.delta_mu_pm_minus_maxwell_MeV, -0.03; atol=1e-10)

    cmp_pm_only = Models._pm_compare_with_maxwell(291.16, nothing; comparison_mu_tol=0.05)
    @test cmp_pm_only.comparison_status == :pm_only

    @test isdefined(Models, :_pm_pressure_crosscheck)
    crosscheck = Models._pm_pressure_crosscheck(rows_transition)
    @test crosscheck.consistent
end

@testset "PM phase diagnostic acceptance and continuity" begin
    @test isdefined(Models, :_pm_accept_branch_point)
    accepted = Models._pm_accept_branch_point((
        converged=true,
        omega=-1.0,
        pressure=1.0,
        rho_norm=0.5,
        residual_norm=1e-7,
    ); residual_accept_tol=1e-6)
    @test accepted.accepted
    @test accepted.branch_status == :accepted

    rejected = Models._pm_accept_branch_point((
        converged=true,
        omega=-1.0,
        pressure=1.0,
        rho_norm=0.5,
        residual_norm=1e-3,
    ); residual_accept_tol=1e-6)
    @test !rejected.accepted
    @test rejected.branch_status == :nonconverged

    invalid_thermo = Models._pm_accept_branch_point((
        converged=true,
        omega=NaN,
        pressure=1.0,
        rho_norm=0.5,
        residual_norm=1e-7,
    ); residual_accept_tol=1e-6)
    @test !invalid_thermo.accepted
    @test invalid_thermo.branch_status == :invalid_thermo

    @test isdefined(Models, :_pm_check_branch_continuity)
    continuity_ok = Models._pm_check_branch_continuity(
        [1.0, 1.0, 1.0],
        [1.1, 1.0, 1.0],
        0.50,
        0.58;
        continuity_x_tol=0.25,
        continuity_rho_tol=0.15,
    )
    @test continuity_ok.continuity_ok
    @test continuity_ok.branch_status == :accepted

    continuity_jump = Models._pm_check_branch_continuity(
        [1.0, 1.0, 1.0],
        [1.5, 1.0, 1.0],
        0.50,
        0.80;
        continuity_x_tol=0.25,
        continuity_rho_tol=0.15,
    )
    @test !continuity_jump.continuity_ok
    @test continuity_jump.branch_status == :branch_jump
    @test continuity_jump.endpoint_cause == :branch_jump
end

@testset "PM phase diagnostic artifact schema helpers" begin
    @test isdefined(Models, :_pm_branch_scan_fieldnames)
    @test isdefined(Models, :_pm_branch_scan_record)

    fieldnames = Models._pm_branch_scan_fieldnames()
    @test fieldnames == (
        :T_MeV,
        :mu_MeV,
        :branch,
        :branch_status,
        :seed_source,
        :continuity_ok,
        :converged,
        :residual_norm,
        :rho_norm,
        :pressure,
        :omega,
        :endpoint_cause,
    )

    branch_record = Models._pm_branch_scan_record((
        T_MeV=130.9,
        mu_MeV=291.0,
        branch=:hadron,
        branch_status=:accepted,
        seed_source=:seed0,
        continuity_ok=true,
        converged=true,
        residual_norm=1e-7,
        rho_norm=0.5,
        pressure=0.12,
        omega=-0.12,
        endpoint_cause=nothing,
    ))
    @test keys(branch_record) == fieldnames
    @test branch_record.endpoint_cause === nothing
end

@testset "PM phase diagnostic Maxwell reference wiring" begin
    @test isdefined(Models, :_pm_maxwell_reference_from_rows)

    branch_rows = [
        (T_MeV=130.9, mu_MeV=290.90, branch=:hadron, accepted=true, rho_norm=0.10, xi=0.0, solver_backend=:models, p_num=24, t_num=8),
        (T_MeV=130.9, mu_MeV=291.00, branch=:hadron, accepted=true, rho_norm=0.20, xi=0.0, solver_backend=:models, p_num=24, t_num=8),
        (T_MeV=130.9, mu_MeV=291.10, branch=:hadron, accepted=true, rho_norm=0.30, xi=0.0, solver_backend=:models, p_num=24, t_num=8),
        (T_MeV=130.9, mu_MeV=291.20, branch=:hadron, accepted=true, rho_norm=0.40, xi=0.0, solver_backend=:models, p_num=24, t_num=8),
        (T_MeV=130.9, mu_MeV=291.22, branch=:hadron, accepted=true, rho_norm=0.50, xi=0.0, solver_backend=:models, p_num=24, t_num=8),
        (T_MeV=130.9, mu_MeV=291.18, branch=:hadron, accepted=true, rho_norm=0.60, xi=0.0, solver_backend=:models, p_num=24, t_num=8),
        (T_MeV=130.9, mu_MeV=291.10, branch=:quark, accepted=true, rho_norm=0.70, xi=0.0, solver_backend=:models, p_num=24, t_num=8),
        (T_MeV=130.9, mu_MeV=291.05, branch=:quark, accepted=true, rho_norm=0.80, xi=0.0, solver_backend=:models, p_num=24, t_num=8),
        (T_MeV=130.9, mu_MeV=291.09, branch=:quark, accepted=true, rho_norm=0.90, xi=0.0, solver_backend=:models, p_num=24, t_num=8),
        (T_MeV=130.9, mu_MeV=291.16, branch=:quark, accepted=true, rho_norm=1.00, xi=0.0, solver_backend=:models, p_num=24, t_num=8),
        (T_MeV=130.9, mu_MeV=291.24, branch=:quark, accepted=true, rho_norm=1.10, xi=0.0, solver_backend=:models, p_num=24, t_num=8),
        (T_MeV=130.9, mu_MeV=291.33, branch=:quark, accepted=true, rho_norm=1.20, xi=0.0, solver_backend=:models, p_num=24, t_num=8),
        (T_MeV=130.9, mu_MeV=291.42, branch=:quark, accepted=true, rho_norm=1.30, xi=0.0, solver_backend=:models, p_num=24, t_num=8),
    ]

    maxwell = Models._pm_maxwell_reference_from_rows(
        branch_rows;
        T_MeV=130.9,
        xi=0.0,
        solver_backend=:models,
        p_num=24,
        t_num=8,
    )
    @test maxwell.mu_transition_maxwell_MeV !== nothing
    @test maxwell.comparison_status == :maxwell_only

    mismatched = Models._pm_maxwell_reference_from_rows(
        branch_rows;
        T_MeV=130.9,
        xi=0.1,
        solver_backend=:models,
        p_num=24,
        t_num=8,
    )
    @test mismatched.mu_transition_maxwell_MeV === nothing
    @test mismatched.comparison_status == :neither
end

@testset "PM phase diagnostic local mu refinement" begin
    @test isdefined(Models, :_pm_refine_transition_bracket)

    coarse_rows = [
        (mu_MeV=291.10, delta_omega=0.02, delta_pressure=-0.02, hadron_exists=true, quark_exists=true),
        (mu_MeV=291.20, delta_omega=-0.01, delta_pressure=0.01, hadron_exists=true, quark_exists=true),
    ]
    refined = Models._pm_refine_transition_bracket(coarse_rows; mu_refine_step=0.01)

    @test length(refined) == 11
    @test refined[1].mu_MeV == 291.10
    @test refined[end].mu_MeV == 291.20
    @test refined[2].mu_MeV == 291.11
end

@testset "PM phase diagnostic input contract" begin
    @test_throws ArgumentError Models.derive_pm_seed_pair(130.9, [291.2, 291.0, 291.1])
    @test_throws ArgumentError Models.analyze_pm_branch_competition(
        T_values=130.9,
        mu_grid=[290.9, 291.0],
        xi=0.0,
        solver_backend=:models,
        p_num=24,
        t_num=8,
        output_dir=tempname(),
    )

    @test_throws ArgumentError Models.analyze_pm_branch_competition(
        T_values=[130.9],
        mu_grid=[290.9, NaN],
        xi=0.0,
        solver_backend=:models,
        p_num=24,
        t_num=8,
        output_dir=tempname(),
    )
end

@testset "PM helper exports stay internal" begin
    @test !Base.isexported(Models, :_pm_interpolate_transition_mu)
    @test !Base.isexported(Models, :_pm_branch_scan_fieldnames)
    @test !Base.isexported(Models, :_pm_maxwell_reference_from_rows)
    @test !Base.isexported(Models, :_pm_refine_transition_bracket)
end
