"""
Sparse diagnostic scan for the preferred candidate route:

    BQS quark-only equilibrium -> legacy/full charged-KMT BU A/B
    -> K+/pi+ and K-/pi- along the baseline chemical freeze-out curve.

This is intentionally separate from the partial-feedback scan.  The meson
densities are post-processing observables and do not enter the BQS stationarity
equations or the PNJL grand potential in this script.
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(@__DIR__, "meson_conserved_charge_feedback_utils.jl"))
include(joinpath(@__DIR__, "meson_conserved_charge_feedback_runtime.jl"))

using JSON
using Printf: @printf
using .Constants_PNJL: ħc_MeV_fm
using .Models
using .MesonConservedChargeFeedbackUtils: charged_meson_conserved_densities,
                                           choose_freezeout_sqrts_grid,
                                           total_conserved_charge_residual
using .MesonConservedChargeFeedbackRuntime: FeedbackSettings, solve_quark_only_bu_ab_point

const FULL_KMT_COST_FACTOR = 2.0

@inline _env_float(name::AbstractString, default::Real) = parse(Float64, get(ENV, name, string(default)))
@inline _env_int(name::AbstractString, default::Integer) = parse(Int, get(ENV, name, string(default)))
@inline _ratio(numerator::Real, denominator::Real) = denominator > 0.0 ? numerator / denominator : NaN

function _settings()
    return FeedbackSettings(
        label="full_kmt_bu_freezeout_diagnostic",
        qmax=_env_float("MESON_FEEDBACK_QMAX", 4.0),
        q_nodes=_env_int("MESON_FEEDBACK_Q_NODES", 4),
        omega_min=_env_float("MESON_FEEDBACK_OMEGA_MIN", 0.05),
        omega_max=_env_float("MESON_FEEDBACK_OMEGA_MAX", 3.0),
        omega_nodes=_env_int("MESON_FEEDBACK_OMEGA_NODES", 8),
        eta=_env_float("MESON_FEEDBACK_ETA", 1e-6),
        density_policy=:x_min_cut,
        bose_x_min=_env_float("MESON_FEEDBACK_BOSE_X_MIN", 0.05),
    )
end

function _benchmark_median_seconds()
    if haskey(ENV, "MESON_FEEDBACK_BENCH_MEDIAN_S")
        value = parse(Float64, ENV["MESON_FEEDBACK_BENCH_MEDIAN_S"])
        isfinite(value) && value > 0.0 || throw(ArgumentError("MESON_FEEDBACK_BENCH_MEDIAN_S must be positive"))
        return value, :environment
    end
    summary_path = get(ENV, "MESON_FEEDBACK_BENCH_SUMMARY", joinpath(
        PROJECT_ROOT,
        "data", "outputs", "perf", "relaxtime", "meson_conserved_charge_partial_feedback",
        "hot_start_summary.json",
    ))
    isfile(summary_path) || throw(ArgumentError(
        "benchmark summary not found at $(summary_path); run the benchmark or set MESON_FEEDBACK_BENCH_MEDIAN_S",
    ))
    parsed = JSON.parsefile(summary_path)
    median_s = Float64(parsed["elapsed_s"]["median"])
    isfinite(median_s) && median_s > 0.0 || throw(ArgumentError("benchmark summary median must be positive"))
    return median_s, Symbol(summary_path)
end

function _grid_from_environment(median_s::Float64)
    if haskey(ENV, "MESON_FEEDBACK_SCAN_SQRTS")
        values = Float64[parse(Float64, strip(x)) for x in split(ENV["MESON_FEEDBACK_SCAN_SQRTS"], ',') if !isempty(strip(x))]
        isempty(values) && throw(ArgumentError("MESON_FEEDBACK_SCAN_SQRTS must not be empty"))
        return sort(unique(values)), :explicit
    end
    budget = _env_float("MESON_FEEDBACK_SCAN_BUDGET_S", 600.0)
    # A/B evaluates the BU kernel twice; reserve a conservative factor over
    # the previous partial-feedback median while keeping the same grid policy.
    return choose_freezeout_sqrts_grid(FULL_KMT_COST_FACTOR * median_s; budget_seconds=budget), :benchmark_budget
end

@inline function _status_string(result)
    return String(result.status)
end

function _meson_payload(density_set)
    return (
        n_pi_plus=Float64(density_set.n_pi_plus),
        n_pi_minus=Float64(density_set.n_pi_minus),
        n_K_plus=Float64(density_set.n_K_plus),
        n_K_minus=Float64(density_set.n_K_minus),
        Kplus_over_piplus=Float64(density_set.Kplus_over_piplus),
        Kminus_over_piminus=Float64(density_set.Kminus_over_piminus),
        pi_plus_status=_status_string(density_set.pi_plus),
        pi_minus_status=_status_string(density_set.pi_minus),
        K_plus_status=_status_string(density_set.K_plus),
        K_minus_status=_status_string(density_set.K_minus),
        pi_plus_min_E_minus_mu=Float64(density_set.pi_plus.min_E_minus_mu),
        pi_minus_min_E_minus_mu=Float64(density_set.pi_minus.min_E_minus_mu),
        K_plus_min_E_minus_mu=Float64(density_set.K_plus.min_E_minus_mu),
        K_minus_min_E_minus_mu=Float64(density_set.K_minus.min_E_minus_mu),
    )
end

function _charge_diagnostic(payload)
    try
        return (value=charged_meson_conserved_densities(payload), message="")
    catch err
        return (
            value=(rho_B=NaN, rho_Q=NaN, rho_S=NaN),
            message="invalid meson density for charge diagnostic: $(sprint(showerror, err))",
        )
    end
end

function _row(
    point,
    output,
    seed_source::Symbol,
    target_ratio::Float64,
    rho_S_target::Float64,
    rho0::Float64,
    settings,
    p_num::Int,
    t_num::Int,
    quark_residual_norm_max::Float64,
)
    if !output.baseline.converged
        return (
            sqrt_s_NN_GeV=point.sqrt_s_NN_GeV,
            T_MeV=point.T_MeV,
            muB_MeV=point.muB_MeV,
            path_profile="baseline_freezeout",
            route="full_kmt_charged_coupling_quark_only",
            production_candidate_status="not_authorized",
            status="quark_only_failed",
            message=String(output.message),
            baseline_converged=false,
            bu_status="not_evaluated",
            quark_only_muQ_MeV=NaN,
            quark_only_muS_MeV=NaN,
            mu_u_MeV=NaN,
            mu_d_MeV=NaN,
            mu_s_MeV=NaN,
            m_u_inv_fm=NaN,
            m_d_inv_fm=NaN,
            m_s_inv_fm=NaN,
            rho_B_q=NaN,
            rho_Q_q=NaN,
            rho_S_q=NaN,
            bqs_residual_norm=NaN,
            rho_Q_M_legacy=NaN,
            rho_S_M_legacy=NaN,
            rho_Q_M_full=NaN,
            rho_S_M_full=NaN,
            rho_Q_total_legacy=NaN,
            rho_S_total_legacy=NaN,
            rho_Q_total_full=NaN,
            rho_S_total_full=NaN,
            residual_norm_legacy=NaN,
            residual_norm_full=NaN,
            n_pi_plus_legacy=NaN,
            n_pi_minus_legacy=NaN,
            n_K_plus_legacy=NaN,
            n_K_minus_legacy=NaN,
            Kplus_over_piplus_legacy=NaN,
            Kminus_over_piminus_legacy=NaN,
            n_pi_plus_full=NaN,
            n_pi_minus_full=NaN,
            n_K_plus_full=NaN,
            n_K_minus_full=NaN,
            Kplus_over_piplus_full=NaN,
            Kminus_over_piminus_full=NaN,
            pi_plus_status_legacy="not_evaluated",
            pi_minus_status_legacy="not_evaluated",
            K_plus_status_legacy="not_evaluated",
            K_minus_status_legacy="not_evaluated",
            pi_plus_status_full="not_evaluated",
            pi_minus_status_full="not_evaluated",
            K_plus_status_full="not_evaluated",
            K_minus_status_full="not_evaluated",
            pi_plus_min_E_minus_mu_full=NaN,
            pi_minus_min_E_minus_mu_full=NaN,
            K_plus_min_E_minus_mu_full=NaN,
            K_minus_min_E_minus_mu_full=NaN,
            K12_P=NaN,
            K45_P=NaN,
            K67_P=NaN,
            K03_P=NaN,
            K38_P=NaN,
            KMT_G_fm2=NaN,
            KMT_K_fm5=NaN,
            seed_source=String(seed_source),
            gap_residual_norm=NaN,
            baseline_elapsed_s=output.baseline_elapsed_s,
            legacy_density_elapsed_s=NaN,
            full_density_elapsed_s=NaN,
            total_elapsed_s=output.total_elapsed_s,
            target_Q_over_B=target_ratio,
            rho_S_target=rho_S_target,
            rho0_fm3=rho0,
            p_num=p_num,
            t_num=t_num,
            quark_residual_norm_max=quark_residual_norm_max,
            qmax_fm_inv=settings.qmax,
            q_nodes=settings.q_nodes,
            omega_min_fm_inv=settings.omega_min,
            omega_max_fm_inv=settings.omega_max,
            omega_nodes=settings.omega_nodes,
            eta_fm_inv=settings.eta,
            density_policy=String(settings.density_policy),
            bose_x_min=settings.bose_x_min,
            scheme="current",
            phase_convention="arg_propagator",
            real_axis_mode="finite_eta",
        )
    end

    baseline_mu = output.baseline_mu
    bqs = output.baseline_bqs
    legacy = _meson_payload(output.legacy)
    full = _meson_payload(output.full)
    legacy_charge_result = _charge_diagnostic(output.legacy)
    full_charge_result = _charge_diagnostic(output.full)
    legacy_charge = legacy_charge_result.value
    full_charge = full_charge_result.value
    target_q = target_ratio * bqs.rho_B
    bqs_residual = sqrt((bqs.rho_Q - target_q)^2 + (bqs.rho_S - rho_S_target)^2)
    total_legacy = isempty(legacy_charge_result.message) ? total_conserved_charge_residual(
        (rho_B=bqs.rho_B, rho_Q=bqs.rho_Q, rho_S=bqs.rho_S),
        legacy_charge;
        charge_to_baryon_ratio=target_ratio,
        strangeness_density_target=rho_S_target,
        rho0=rho0,
    ) : (rho_Q_total=NaN, rho_S_total=NaN, norm=NaN)
    total_full = isempty(full_charge_result.message) ? total_conserved_charge_residual(
        (rho_B=bqs.rho_B, rho_Q=bqs.rho_Q, rho_S=bqs.rho_S),
        full_charge;
        charge_to_baryon_ratio=target_ratio,
        strangeness_density_target=rho_S_target,
        rho0=rho0,
    ) : (rho_Q_total=NaN, rho_S_total=NaN, norm=NaN)
    kc = output.kernel_couplings
    invalid_charge = !isempty(legacy_charge_result.message) || !isempty(full_charge_result.message)
    status = output.converged && !invalid_charge ? "quark_only_full_kmt_bu_ok" : "bu_density_status_not_ok"
    message = isempty(legacy_charge_result.message) && isempty(full_charge_result.message) ?
        String(output.message) : join(filter(!isempty, (legacy_charge_result.message, full_charge_result.message)), "; ")
    return (
        sqrt_s_NN_GeV=point.sqrt_s_NN_GeV,
        T_MeV=point.T_MeV,
        muB_MeV=point.muB_MeV,
        path_profile="baseline_freezeout",
        route="full_kmt_charged_coupling_quark_only",
        production_candidate_status="not_authorized",
        status=status,
        message=message,
        baseline_converged=true,
        bu_status=output.converged ? "ok" : String(output.reason),
        quark_only_muQ_MeV=baseline_mu.mu_Q * ħc_MeV_fm,
        quark_only_muS_MeV=baseline_mu.mu_S * ħc_MeV_fm,
        mu_u_MeV=output.qp.μ.u * ħc_MeV_fm,
        mu_d_MeV=output.qp.μ.d * ħc_MeV_fm,
        mu_s_MeV=output.qp.μ.s * ħc_MeV_fm,
        m_u_inv_fm=output.masses.u,
        m_d_inv_fm=output.masses.d,
        m_s_inv_fm=output.masses.s,
        rho_B_q=bqs.rho_B,
        rho_Q_q=bqs.rho_Q,
        rho_S_q=bqs.rho_S,
        bqs_residual_norm=bqs_residual,
        rho_Q_M_legacy=legacy_charge.rho_Q,
        rho_S_M_legacy=legacy_charge.rho_S,
        rho_Q_M_full=full_charge.rho_Q,
        rho_S_M_full=full_charge.rho_S,
        rho_Q_total_legacy=total_legacy.rho_Q_total,
        rho_S_total_legacy=total_legacy.rho_S_total,
        rho_Q_total_full=total_full.rho_Q_total,
        rho_S_total_full=total_full.rho_S_total,
        residual_norm_legacy=total_legacy.norm,
        residual_norm_full=total_full.norm,
        n_pi_plus_legacy=legacy.n_pi_plus,
        n_pi_minus_legacy=legacy.n_pi_minus,
        n_K_plus_legacy=legacy.n_K_plus,
        n_K_minus_legacy=legacy.n_K_minus,
        Kplus_over_piplus_legacy=legacy.Kplus_over_piplus,
        Kminus_over_piminus_legacy=legacy.Kminus_over_piminus,
        n_pi_plus_full=full.n_pi_plus,
        n_pi_minus_full=full.n_pi_minus,
        n_K_plus_full=full.n_K_plus,
        n_K_minus_full=full.n_K_minus,
        Kplus_over_piplus_full=full.Kplus_over_piplus,
        Kminus_over_piminus_full=full.Kminus_over_piminus,
        pi_plus_status_legacy=legacy.pi_plus_status,
        pi_minus_status_legacy=legacy.pi_minus_status,
        K_plus_status_legacy=legacy.K_plus_status,
        K_minus_status_legacy=legacy.K_minus_status,
        pi_plus_status_full=full.pi_plus_status,
        pi_minus_status_full=full.pi_minus_status,
        K_plus_status_full=full.K_plus_status,
        K_minus_status_full=full.K_minus_status,
        pi_plus_min_E_minus_mu_full=full.pi_plus_min_E_minus_mu,
        pi_minus_min_E_minus_mu_full=full.pi_minus_min_E_minus_mu,
        K_plus_min_E_minus_mu_full=full.K_plus_min_E_minus_mu,
        K_minus_min_E_minus_mu_full=full.K_minus_min_E_minus_mu,
        K12_P=kc.K12_P,
        K45_P=kc.K45_P,
        K67_P=kc.K67_P,
        K03_P=kc.K03_P,
        K38_P=kc.K38_P,
        KMT_G_fm2=kc.G,
        KMT_K_fm5=kc.K,
        seed_source=String(seed_source),
        gap_residual_norm=Float64(output.baseline.residual_norm),
        baseline_elapsed_s=output.baseline_elapsed_s,
        legacy_density_elapsed_s=output.legacy_density_elapsed_s,
        full_density_elapsed_s=output.full_density_elapsed_s,
        total_elapsed_s=output.total_elapsed_s,
        target_Q_over_B=target_ratio,
        rho_S_target=rho_S_target,
        rho0_fm3=rho0,
        p_num=p_num,
        t_num=t_num,
        quark_residual_norm_max=quark_residual_norm_max,
        qmax_fm_inv=settings.qmax,
        q_nodes=settings.q_nodes,
        omega_min_fm_inv=settings.omega_min,
        omega_max_fm_inv=settings.omega_max,
        omega_nodes=settings.omega_nodes,
        eta_fm_inv=settings.eta,
        density_policy=String(settings.density_policy),
        bose_x_min=settings.bose_x_min,
        scheme="current",
        phase_convention="arg_propagator",
        real_axis_mode="finite_eta",
    )
end

function _write_csv(path::String, rows)
    mkpath(dirname(path))
    columns = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(String.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, col)) for col in columns), ','))
        end
    end
end

function _write_summary(path::String, rows, median_s::Float64, grid_source, grid, settings)
    mkpath(dirname(path))
    successful = filter(row -> row.status == "quark_only_full_kmt_bu_ok", rows)
    open(path, "w") do io
        println(io, "# Full KMT charged BU A/B freeze-out sparse scan")
        println(io)
        println(io, "- status: diagnostic-only; production candidate authorization: `not_authorized`")
        println(io, "- route: full charged `K12/K45` coupling on the quark-only BQS background")
        println(io, "- freeze-out profile: `default` / `baseline_freezeout`")
        println(io, "- benchmark median reference: `$(median_s) s`; A/B cost factor: `$(FULL_KMT_COST_FACTOR)`")
        println(io, "- grid source: `$(grid_source)`")
        println(io, "- grid (GeV): `$(join(grid, ", "))`")
        println(io, "- density policy: `$(settings.density_policy)`, `bose_x_min=$(settings.bose_x_min)`")
        println(io, "- no `Omega_M` feedback, strict-support or node-convergence claim is made")
        println(io)
        println(io, "| sqrt(s_NN) [GeV] | T [MeV] | status | legacy K+/pi+ | full K+/pi+ | legacy K-/pi- | full K-/pi- | full residual | time [s] |")
        println(io, "|---:|---:|---|---:|---:|---:|---:|---:|---:|")
        for row in sort(rows; by=r -> r.sqrt_s_NN_GeV)
            println(io, "| $(row.sqrt_s_NN_GeV) | $(row.T_MeV) | $(row.status) | $(row.Kplus_over_piplus_legacy) | $(row.Kplus_over_piplus_full) | $(row.Kminus_over_piminus_legacy) | $(row.Kminus_over_piminus_full) | $(row.residual_norm_full) | $(row.total_elapsed_s) |")
        end
        println(io)
        println(io, "successful_points=$(length(successful))/$(length(rows))")
        if !isempty(successful)
            plus = [row.Kplus_over_piplus_full for row in successful if isfinite(row.Kplus_over_piplus_full)]
            minus = [row.Kminus_over_piminus_full for row in successful if isfinite(row.Kminus_over_piminus_full)]
            if !isempty(plus) && !isempty(minus)
                println(io, "full_ratio_range_Kplus_over_piplus=$(minimum(plus))..$(maximum(plus))")
                println(io, "full_ratio_range_Kminus_over_piminus=$(minimum(minus))..$(maximum(minus))")
            end
        end
        println(io)
        println(io, "本文件只支持同一 quark-only 背景下的 legacy/full charged KMT BU A/B 观察；任何生产晋升仍需 strict-support、节点收敛、外部参考和作者审查。")
    end
end

function main()
    median_s, median_source = _benchmark_median_seconds()
    grid, grid_source = _grid_from_environment(median_s)
    settings = _settings()
    target_ratio = _env_float("MESON_FEEDBACK_TARGET_QB", 0.4)
    rho_S_target = _env_float("MESON_FEEDBACK_RHOS_TARGET", 0.0)
    rho0 = _env_float("MESON_FEEDBACK_RHO0", 0.16)
    profile = Models.load_freezeout_profile(profile="default")
    points = Models.build_freezeout_scan_points(grid; profile=profile, traversal=:sqrts_descending)

    rows = NamedTuple[]
    previous_solution = nothing
    p_num = _env_int("MESON_FEEDBACK_P_NUM", 8)
    t_num = _env_int("MESON_FEEDBACK_T_NUM", 4)
    quark_residual_norm_max = _env_float("MESON_FEEDBACK_QUARK_RESIDUAL_MAX", 1e-5)
    for point in points
        seed_source = previous_solution === nothing ? :cold_or_internal : :previous_quark_only_solution
        model = create_model(:PNJL)
        output = solve_quark_only_bu_ab_point(
            model,
            point.T_fm,
            point.muB_fm,
            settings;
            baseline_seed=previous_solution,
            target_ratio=target_ratio,
            rho_S_target=rho_S_target,
            rho0=rho0,
            p_num=p_num,
            t_num=t_num,
            quark_residual_norm_max=quark_residual_norm_max,
            quark_iterations=_env_int("MESON_FEEDBACK_QUARK_ITERATIONS", 200),
        )
        row = _row(
            point,
            output,
            seed_source,
            target_ratio,
            rho_S_target,
            rho0,
            settings,
            p_num,
            t_num,
            quark_residual_norm_max,
        )
        push!(rows, row)
        if output.baseline.converged
            previous_solution = output.baseline.solution
        end
        @printf("[freezeout-full-kmt-ab] sqrt(s)=%.3f GeV status=%s time=%.3fs\n",
                point.sqrt_s_NN_GeV, row.status, row.total_elapsed_s)
    end

    output_dir = get(ENV, "MESON_FEEDBACK_SCAN_OUTPUT_DIR", joinpath(
        PROJECT_ROOT,
        "data", "outputs", "results", "relaxtime", "analysis",
        "full_kmt_bu_freezeout_sparse",
    ))
    csv_path = joinpath(output_dir, "full_kmt_bu_ab_freezeout_scan.csv")
    summary_path = joinpath(output_dir, "README.md")
    _write_csv(csv_path, sort(rows; by=r -> r.sqrt_s_NN_GeV))
    _write_summary(summary_path, rows, median_s, median_source, grid, settings)
    println("[freezeout-full-kmt-ab] csv=$(csv_path)")
    println("[freezeout-full-kmt-ab] summary=$(summary_path)")
    return (rows=rows, csv=csv_path, summary=summary_path, grid=grid)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
