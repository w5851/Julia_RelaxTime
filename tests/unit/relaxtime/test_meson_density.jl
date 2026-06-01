"""
MesonDensity 模块单元测试（稳定粒子极限）。
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

const _RELAXTIME_PATH_MD = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAXTIME_PATH_MD)
end

using Main.RelaxTime.MesonDensity: DEFAULT_MESON_DENSITY_Q_NODES,
                                   DEFAULT_PHASE_SHIFT_OMEGA_NODES,
                                   DEFAULT_PHASE_SHIFT_OMEGA_MAX,
                                   DEFAULT_PHASE_SHIFT_Q_MAX,
                                   DEFAULT_PHASE_SHIFT_Q_NODES,
                                   _phase_shift_scheme_symbol,
                                   _real_axis_mode_symbol,
                                   _density_policy_symbol,
                                   _phase_display_symbol,
                                   _apply_phase_display,
                                   _noanom_policy_symbol,
                                   _apply_noanom_policy,
                                   _unwrap_phases,
                                   bose_distribution,
                                   meson_degeneracy,
                                   phase_shift_point_diagnostic,
                                   phase_shift_meson_density_derivative_reference_summary,
                                   phase_shift_meson_number_density_derivative_reference,
                                   phase_shift_meson_number_density,
                                   phase_shift_meson_density_summary,
                                   strict_bw_meson_density_summary,
                                   strict_bw_meson_number_density,
                                   stable_meson_number_density,
                                   stable_kpi_ratio,
                                   stable_kpi_scan

const ZETA3 = 1.2020569031595942

@testset "MesonDensity 简并因子与玻色分布" begin
    @test meson_degeneracy(:pi) == 3
    @test meson_degeneracy(:K) == 4
    @test meson_degeneracy(:pi_plus) == 1
    @test meson_degeneracy(:pi_minus) == 1
    @test meson_degeneracy(:K_plus) == 1
    @test meson_degeneracy(:K_minus) == 1
    @test meson_degeneracy(:sigma_pi) == 1
    @test meson_degeneracy(:sigma_K) == 4
    @test meson_degeneracy(:pi; charge_resolved=true) == 1
    @test meson_degeneracy(:K; charge_resolved=true) == 1

    @test bose_distribution(1.0, 0.0, 0.2) > 0.0
    @test bose_distribution(1.0, 0.0, 0.0) == 0.0
    @test_throws ArgumentError bose_distribution(0.5, 0.5, 0.2)
end

@testset "MesonDensity 稳定粒子极限" begin
    T = 0.25
    n_pi = stable_meson_number_density(0.14, T; degeneracy=3, num_q_nodes=DEFAULT_MESON_DENSITY_Q_NODES)
    n_K = stable_meson_number_density(0.49, T; degeneracy=4, num_q_nodes=DEFAULT_MESON_DENSITY_Q_NODES)

    @test n_pi > 0.0
    @test n_K > 0.0
    @test n_pi > n_K

    n_low = stable_meson_number_density(0.14, 0.15; degeneracy=3)
    n_high = stable_meson_number_density(0.14, 0.20; degeneracy=3)
    @test n_high > n_low

    # 质量为零、μ=0 时的解析结果：n = d * ζ(3) / π² * T³
    n_massless = stable_meson_number_density(0.0, 0.3; degeneracy=2, μ=-1e-6, qmax=12.0, num_q_nodes=512)
    n_expected = 2.0 * ZETA3 / π^2 * 0.3^3
    @test n_massless ≈ n_expected rtol=2e-3

    @test_throws ArgumentError stable_meson_number_density(0.14, 0.2; μ=0.14)
end

@testset "MesonDensity Phase-E3 参数校验" begin
    @test DEFAULT_PHASE_SHIFT_Q_MAX == 12.0
    @test DEFAULT_PHASE_SHIFT_Q_NODES == 48
    @test DEFAULT_PHASE_SHIFT_OMEGA_MAX == 10.0
    @test DEFAULT_PHASE_SHIFT_OMEGA_NODES == 48
    @test _phase_shift_scheme_symbol(:current) == :phase_shift_current
    @test _phase_shift_scheme_symbol(:gbu) == :phase_shift_gbu_reference
    @test _real_axis_mode_symbol(:finite_eta) == :finite_eta
    @test _real_axis_mode_symbol(:bu2020_pv_eta0) == :pv_b0_eta0
    @test _density_policy_symbol(:strict) == :strict_normal_domain
    @test _density_policy_symbol(:excitation_only) == :excitation_only_E_gt_mu
    @test _phase_display_symbol(:unwrapped) == :unwrapped
    @test _phase_display_symbol(:folded_0_pi) == :fold_0_pi
    @test _noanom_policy_symbol(:none) == :none
    @test _noanom_policy_symbol(:temp7_low_energy_branch_subtraction) == :low_energy_branch_subtraction

    qp = (m=(u=1.0, d=1.0, s=1.2), μ=(u=0.0, d=0.0, s=0.0), A=(u=0.1, d=0.1, s=0.1))
    tp_bad = (T=0.2, Φ=0.5, Φbar=0.5, ξ=0.1)
    @test_throws ArgumentError phase_shift_meson_density_summary(qp, tp_bad)
end

@testset "MesonDensity phase display policy" begin
    phases = [-π, -0.1, 0.2, 4π + 0.1]
    unwrapped = _apply_phase_display(phases, :unwrapped)
    folded = _apply_phase_display(phases, :fold_0_pi)
    @test unwrapped == Float64.(phases)
    @test all(0.0 .<= folded .<= π)
    @test folded[1] ≈ π
    @test folded[2] ≈ 0.1
    @test folded[4] ≈ 0.1
    @test_throws ArgumentError _apply_phase_display(phases, :principal)
end

@testset "MesonDensity unwrap phases defaults to strict pi threshold" begin
    phases = [0.0, π - 2e-5, π - 1e-5]
    unwrapped = _unwrap_phases(phases)
    tolerant = _unwrap_phases(phases; branch_tol=1e-4)
    @test unwrapped[2] > 3.0
    @test unwrapped[3] > unwrapped[2]
    @test tolerant[2] < -3.0
    @test tolerant[3] > tolerant[2]
end

@testset "MesonDensity unwrap phases rejects negative tolerance" begin
    phases = [0.1, 0.2]
    @test_throws ArgumentError _unwrap_phases(phases; branch_tol=-1e-6)
end

@testset "MesonDensity K/π 扫描" begin
    ratio = stable_kpi_ratio(0.14, 0.49, 0.20)
    @test 0.0 < ratio < 1.0

    scan = stable_kpi_scan([0.16, 0.18, 0.20]; m_pi=0.14, m_K=0.49)
    @test length(scan.temperatures) == 3
    @test length(scan.n_pi) == 3
    @test length(scan.n_K) == 3
    @test length(scan.kpi_ratio) == 3
    @test all(diff(scan.n_pi) .> 0.0)
    @test all(diff(scan.n_K) .> 0.0)
end

@testset "MesonDensity reduced strict BW" begin
    T = 0.20
    stable = stable_meson_number_density(0.14, T; degeneracy=3, qmax=12.0, num_q_nodes=48)
    bw_zero = strict_bw_meson_number_density(
        0.14,
        0.0,
        T;
        degeneracy=3,
        qmax=12.0,
        q_nodes=48,
        omega_max=10.0,
        omega_nodes=24,
    )
    @test bw_zero.mode == :stable_limit
    @test bw_zero.density ≈ stable rtol=1e-10

    bw_finite = strict_bw_meson_number_density(
        0.14,
        0.08,
        T;
        degeneracy=3,
        qmax=12.0,
        q_nodes=24,
        omega_max=10.0,
        omega_nodes=24,
    )
    @test bw_finite.mode == :strict_bw_reduced
    @test isfinite(bw_finite.density)
    @test bw_finite.density > 0.0
    @test isfinite(bw_finite.omega_shell_at_qmax)

    summary = strict_bw_meson_density_summary(
        0.14, 0.08,
        0.49, 0.12,
        T;
        qmax=12.0,
        q_nodes=24,
        omega_max=10.0,
        omega_nodes=24,
    )
    @test summary.n_pi > 0.0
    @test summary.n_K > 0.0
    @test 0.0 < summary.kpi_ratio < 1.5
end

@testset "MesonDensity generalized BU reference scheme" begin
    qp = (m=(u=0.098, d=0.098, s=0.42), μ=(u=0.0, d=0.0, s=0.0), A=(u=0.1, d=0.1, s=0.08))
    tp = (T=0.18, Φ=0.25, Φbar=0.25, ξ=0.0)

    current = phase_shift_meson_number_density(
        :pi_plus,
        qp,
        tp;
        scheme=:current,
        degeneracy=1,
        qmax=4.0,
        q_nodes=6,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=6,
    )
    gbu = phase_shift_meson_number_density(
        :pi_plus,
        qp,
        tp;
        scheme=:gbu_reference,
        degeneracy=1,
        qmax=4.0,
        q_nodes=6,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=6,
    )

    @test current.scheme == :phase_shift_current
    @test gbu.scheme == :phase_shift_gbu_reference
    @test isfinite(current.density)
    @test isfinite(gbu.density)
    @test current.density > 0.0
    @test gbu.density > 0.0
    @test !isapprox(current.density, gbu.density; rtol=1e-6, atol=1e-10)

    summary = phase_shift_meson_density_summary(
        qp,
        tp;
        scheme=:gbu_reference,
        qmax=4.0,
        q_nodes=6,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=6,
    )
    @test summary.scheme == :phase_shift_gbu_reference
    @test summary.pi_density.scheme == :phase_shift_gbu_reference
    @test summary.k_density.scheme == :phase_shift_gbu_reference
    @test summary.n_pi > 0.0
    @test summary.n_K > 0.0
end

@testset "MesonDensity real-axis mode and Bose-domain policy metadata" begin
    qp = (m=(u=0.098, d=0.098, s=0.42), μ=(u=0.0, d=0.0, s=0.0), A=(u=0.1, d=0.1, s=0.08))
    tp = (T=0.18, Φ=0.25, Φbar=0.25, ξ=0.0)

    @test_throws ArgumentError phase_shift_meson_number_density(
        :pi_plus,
        qp,
        tp;
        eta=0.0,
        real_axis_mode=:finite_eta,
        qmax=4.0,
        q_nodes=4,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=4,
    )

    pv = phase_shift_meson_number_density(
        :pi_plus,
        qp,
        tp;
        real_axis_mode=:pv_b0_eta0,
        phase_convention=:arg_inverse_propagator,
        qmax=4.0,
        q_nodes=4,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=4,
    )
    @test pv.real_axis_mode == :pv_b0_eta0
    @test pv.polarization_backend == :pv_b0_real_axis
    @test pv.phase_convention == :arg_inverse_propagator
    @test pv.eta == 0.0
    @test pv.status == :ok

    strict = phase_shift_meson_number_density(
        :pi_plus,
        qp,
        tp;
        μ=0.5,
        qmax=4.0,
        q_nodes=4,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=4,
        density_policy=:strict_normal_domain,
    )
    @test strict.status == :unsafe_bose_domain
    @test isnan(strict.density)
    @test strict.unsafe_bose_count > 0
    @test strict.min_E_minus_mu < 0.0

    boundary_strict = phase_shift_meson_number_density(
        :pi_plus,
        qp,
        tp;
        μ=0.06,
        qmax=4.0,
        q_nodes=4,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=2,
        density_policy=:strict_normal_domain,
    )
    @test boundary_strict.status == :unsafe_bose_domain
    @test boundary_strict.unsafe_bose_count >= 1

    excitation_only = phase_shift_meson_number_density(
        :pi_plus,
        qp,
        tp;
        μ=0.5,
        qmax=4.0,
        q_nodes=4,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=4,
        density_policy=:excitation_only_E_gt_mu,
    )
    @test excitation_only.status == :ok
    @test excitation_only.omega_min > 0.5
    @test isfinite(excitation_only.density)

    noanom = phase_shift_meson_number_density(
        :K_plus,
        qp,
        tp;
        qmax=4.0,
        q_nodes=4,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=4,
        noanom_policy=:low_energy_branch_subtraction,
    )
    @test noanom.noanom_policy == :low_energy_branch_subtraction
    @test noanom.noanom_applied
    @test isfinite(noanom.noanom_landau_omega_min)
    @test isfinite(noanom.noanom_landau_omega_max)
    @test noanom.noanom_landau_omega_max > noanom.noanom_landau_omega_min

    qp_noanom = (m=(u=1.0, d=1.0, s=1.4), μ=(u=0.2, d=0.2, s=0.0), A=(u=0.1, d=0.1, s=0.1))
    omega_grid = [0.1, 0.2, 0.3, 0.4, 0.7, 0.8, 0.9]
    phase_grid = [0.0, 0.2, 0.3, 0.1, 0.0, 1.0, 1.1]
    phases_noanom, diag = _apply_noanom_policy(
        :K_plus,
        omega_grid,
        phase_grid,
        qp_noanom;
        noanom_policy=:low_energy_branch_subtraction,
    )
    @test diag.noanom_applied
    @test diag.noanom_removed_component_count == 1
    @test phases_noanom[2:4] == zeros(3)
    @test phases_noanom[6] > 0.0
end

@testset "MesonDensity derivative-reference helper is AD-backed" begin
    qp = (m=(u=0.098, d=0.098, s=0.42), μ=(u=0.0, d=0.0, s=0.0), A=(u=0.1, d=0.1, s=0.08))
    tp = (T=0.18, Φ=0.25, Φbar=0.25, ξ=0.0)

    single = phase_shift_meson_number_density_derivative_reference(
        :pi, qp, tp;
        scheme=:current,
        qmax=4.0,
        q_nodes=6,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=6,
    )
    @test single.derivative_backend == :forwarddiff
    @test isfinite(single.density)
    @test single.max_formula_abs_diff < 1e-8

    summary = phase_shift_meson_density_derivative_reference_summary(
        qp, tp;
        scheme=:gbu_reference,
        qmax=4.0,
        q_nodes=6,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=6,
    )
    @test summary.derivative_backend == :forwarddiff
    @test summary.scheme == :phase_shift_gbu_reference
    @test isfinite(summary.n_pi)
    @test isfinite(summary.n_K)
    @test summary.max_formula_abs_diff < 1e-8
end

@testset "MesonDensity phase-shift point diagnostic AD vs FD" begin
    qp = (m=(u=0.098, d=0.098, s=0.42), μ=(u=0.0, d=0.0, s=0.0), A=(u=0.1, d=0.1, s=0.08))
    tp_iso = (T=0.18, Φ=0.25, Φbar=0.25, ξ=0.0)
    tp_aniso = (T=0.18, Φ=0.25, Φbar=0.25, ξ=0.1)

    diag_iso = phase_shift_point_diagnostic(:pi, 0.4, 0.5, qp, tp_iso; scheme=:current, fd_step=1e-5)
    @test isfinite(diag_iso.dphase_ad)
    @test isfinite(diag_iso.dphase_fd)
    @test isfinite(diag_iso.dphase_raw_fd)
    @test isfinite(diag_iso.dweighted_ad)
    @test isfinite(diag_iso.dweighted_fd)
    @test diag_iso.dphase_formula_abs_diff < 1e-8
    @test abs(diag_iso.dReD_domega - diag_iso.dReD_fd) > 1e-3
    @test abs(diag_iso.dImD_domega - diag_iso.dImD_fd) > 1e-3

    diag_aniso = phase_shift_point_diagnostic(:pi, 0.4, 0.5, qp, tp_aniso; scheme=:current, fd_step=1e-5)
    @test diag_aniso.xi ≈ 0.1
    @test isfinite(diag_aniso.dphase_ad)
    @test isfinite(diag_aniso.dphase_fd)
    @test isfinite(diag_aniso.dphase_raw_fd)
    @test diag_aniso.dphase_formula_abs_diff < 1e-8
    @test abs(diag_aniso.dReD_domega - diag_aniso.dReD_fd) > 1e-3
    @test abs(diag_aniso.dImD_domega - diag_aniso.dImD_fd) > 1e-3
end
