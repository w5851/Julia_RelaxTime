"""
Sparse diagnostic scan of charged-meson partial feedback on the chemical
freeze-out baseline.

This scan is deliberately separate from `Models.run_freezeout_meson_density_scan`:
the latter remains the stable fixed-profile workflow, while this script solves
the outer `(mu_Q,mu_S)` correction point by point and labels the result as
partial-feedback diagnostic data.
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
using .MesonConservedChargeFeedbackUtils
using .MesonConservedChargeFeedbackRuntime: FeedbackSettings, solve_partial_feedback_point

@inline _env_float(name::AbstractString, default::Real) = parse(Float64, get(ENV, name, string(default)))
@inline _env_int(name::AbstractString, default::Integer) = parse(Int, get(ENV, name, string(default)))
@inline _ratio(numerator::Real, denominator::Real) = denominator > 0.0 ? numerator / denominator : NaN

function _settings()
    return FeedbackSettings(
        label="freezeout_sparse_diagnostic",
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
        "benchmark summary not found at $(summary_path); run the hot-start benchmark first or set MESON_FEEDBACK_BENCH_MEDIAN_S",
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
    return choose_freezeout_sqrts_grid(median_s; budget_seconds=budget), :benchmark_budget
end

function _row(
    point,
    output,
    outer_initial_mu_source::Symbol,
)
    base = output.baseline
    if !base.converged
        return (
            sqrt_s_NN_GeV=point.sqrt_s_NN_GeV,
            T_MeV=point.T_MeV,
            muB_MeV=point.muB_MeV,
            path_profile="baseline_freezeout",
            status="quark_only_failed",
            message=String(output.message),
            baseline_converged=false,
            feedback_converged=false,
            quark_only_muQ_MeV=NaN,
            quark_only_muS_MeV=NaN,
            feedback_muQ_MeV=NaN,
            feedback_muS_MeV=NaN,
            mu_u_MeV=NaN,
            mu_d_MeV=NaN,
            mu_s_MeV=NaN,
            rho_B_q=NaN,
            rho_Q_q=NaN,
            rho_S_q=NaN,
            rho_B_M=0.0,
            rho_Q_M=NaN,
            rho_S_M=NaN,
            rho_B_total=NaN,
            rho_Q_total=NaN,
            rho_S_total=NaN,
            residual_norm=NaN,
            gap_residual_norm=NaN,
            n_pi_plus=NaN,
            n_pi_minus=NaN,
            n_K_plus=NaN,
            n_K_minus=NaN,
            Kplus_over_piplus=NaN,
            Kminus_over_piminus=NaN,
            pi_plus_status="not_evaluated",
            pi_minus_status="not_evaluated",
            K_plus_status="not_evaluated",
            K_minus_status="not_evaluated",
            pi_plus_min_E_minus_mu=NaN,
            pi_minus_min_E_minus_mu=NaN,
            K_plus_min_E_minus_mu=NaN,
            K_minus_min_E_minus_mu=NaN,
            outer_initial_mu_source=String(outer_initial_mu_source),
            outer_iterations=0,
            outer_evaluations=0,
            unique_candidate_count=0,
            baseline_elapsed_s=output.baseline_elapsed_s,
            outer_elapsed_s=NaN,
            gap_elapsed_s=NaN,
            density_elapsed_s=NaN,
            total_elapsed_s=output.total_elapsed_s,
            density_policy="x_min_cut",
            q_nodes=0,
            omega_nodes=0,
            bose_x_min=NaN,
        )
    end

    feedback = output.feedback
    result = feedback.result
    payload = result.payload
    meson = charged_meson_conserved_densities(payload)
    total = total_conserved_charge_residual(
        (rho_B=payload.rho_B_q, rho_Q=payload.rho_Q_q, rho_S=payload.rho_S_q),
        meson;
        charge_to_baryon_ratio=_env_float("MESON_FEEDBACK_TARGET_QB", 0.4),
        strangeness_density_target=_env_float("MESON_FEEDBACK_RHOS_TARGET", 0.0),
        rho0=_env_float("MESON_FEEDBACK_RHO0", 0.16),
    )
    initial_mu = output.baseline_mu
    return (
        sqrt_s_NN_GeV=point.sqrt_s_NN_GeV,
        T_MeV=point.T_MeV,
        muB_MeV=point.muB_MeV,
        path_profile="baseline_freezeout",
        status=result.converged ? "partial_feedback_converged" : "partial_feedback_failed",
        message=String(result.message),
        baseline_converged=true,
        feedback_converged=result.converged,
        quark_only_muQ_MeV=initial_mu.mu_Q * ħc_MeV_fm,
        quark_only_muS_MeV=initial_mu.mu_S * ħc_MeV_fm,
        feedback_muQ_MeV=result.mu_Q * ħc_MeV_fm,
        feedback_muS_MeV=result.mu_S * ħc_MeV_fm,
        mu_u_MeV=payload.mu_u * ħc_MeV_fm,
        mu_d_MeV=payload.mu_d * ħc_MeV_fm,
        mu_s_MeV=payload.mu_s * ħc_MeV_fm,
        rho_B_q=payload.rho_B_q,
        rho_Q_q=payload.rho_Q_q,
        rho_S_q=payload.rho_S_q,
        rho_B_M=0.0,
        rho_Q_M=meson.rho_Q,
        rho_S_M=meson.rho_S,
        rho_B_total=total.rho_B_total,
        rho_Q_total=total.rho_Q_total,
        rho_S_total=total.rho_S_total,
        residual_norm=total.norm,
        gap_residual_norm=payload.gap_residual_norm,
        n_pi_plus=payload.n_pi_plus,
        n_pi_minus=payload.n_pi_minus,
        n_K_plus=payload.n_K_plus,
        n_K_minus=payload.n_K_minus,
        Kplus_over_piplus=_ratio(payload.n_K_plus, payload.n_pi_plus),
        Kminus_over_piminus=_ratio(payload.n_K_minus, payload.n_pi_minus),
        pi_plus_status=String(payload.pi_plus_status),
        pi_minus_status=String(payload.pi_minus_status),
        K_plus_status=String(payload.K_plus_status),
        K_minus_status=String(payload.K_minus_status),
        pi_plus_min_E_minus_mu=payload.pi_plus_min_E_minus_mu,
        pi_minus_min_E_minus_mu=payload.pi_minus_min_E_minus_mu,
        K_plus_min_E_minus_mu=payload.K_plus_min_E_minus_mu,
        K_minus_min_E_minus_mu=payload.K_minus_min_E_minus_mu,
        outer_initial_mu_source=String(outer_initial_mu_source),
        outer_iterations=result.iterations,
        outer_evaluations=result.evaluation_count,
        unique_candidate_count=feedback.unique_candidate_count,
        baseline_elapsed_s=output.baseline_elapsed_s,
        outer_elapsed_s=feedback.outer_elapsed_s,
        gap_elapsed_s=feedback.timing.gap_elapsed_s,
        density_elapsed_s=feedback.timing.density_elapsed_s,
        total_elapsed_s=output.total_elapsed_s,
        density_policy=String(feedback.settings.density_policy),
        q_nodes=feedback.settings.q_nodes,
        omega_nodes=feedback.settings.omega_nodes,
        bose_x_min=feedback.settings.bose_x_min,
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
    open(path, "w") do io
        println(io, "# Partial-feedback chemical freeze-out sparse scan")
        println(io)
        println(io, "- status: diagnostic-only")
        println(io, "- freeze-out profile: `default` / `baseline_freezeout`")
        println(io, "- benchmark median single-point time: `$(median_s) s`")
        println(io, "- grid source: `$(grid_source)`")
        println(io, "- grid (GeV): `$(join(grid, ", "))`")
        println(io, "- density policy: `$(settings.density_policy)`, `bose_x_min=$(settings.bose_x_min)`")
        println(io)
        println(io, "| sqrt(s_NN) [GeV] | T [MeV] | muB [MeV] | status | K+/pi+ | K-/pi- | residual | time [s] |")
        println(io, "|---:|---:|---:|---|---:|---:|---:|---:|")
        for row in sort(rows; by=r -> r.sqrt_s_NN_GeV)
            println(io, "| $(row.sqrt_s_NN_GeV) | $(row.T_MeV) | $(row.muB_MeV) | $(row.status) | $(row.Kplus_over_piplus) | $(row.Kminus_over_piminus) | $(row.residual_norm) | $(row.total_elapsed_s) |")
        end
        println(io)
        println(io, "本文件只用于 partial-feedback 的趋势与成本诊断，不承担 strict-support、节点收敛、完整热力学反馈或实验逐点拟合结论。")
    end
end

function main()
    median_s, median_source = _benchmark_median_seconds()
    grid, grid_source = _grid_from_environment(median_s)
    settings = _settings()
    profile = Models.load_freezeout_profile(profile="default")
    points = Models.build_freezeout_scan_points(grid; profile=profile, traversal=:sqrts_descending)

    rows = NamedTuple[]
    previous_solution = nothing
    previous_feedback = nothing
    for point in points
        initial_mu_Q = previous_feedback === nothing ? nothing : previous_feedback.mu_Q
        initial_mu_S = previous_feedback === nothing ? nothing : previous_feedback.mu_S
        initial_source = previous_feedback === nothing ? :current_quark_only : :previous_feedback
        model = create_model(:PNJL)
        output = solve_partial_feedback_point(
            model,
            point.T_fm,
            point.muB_fm,
            settings;
            baseline_seed=previous_solution,
            initial_mu_Q=initial_mu_Q,
            initial_mu_S=initial_mu_S,
            target_ratio=_env_float("MESON_FEEDBACK_TARGET_QB", 0.4),
            rho_S_target=_env_float("MESON_FEEDBACK_RHOS_TARGET", 0.0),
            rho0=_env_float("MESON_FEEDBACK_RHO0", 0.16),
            p_num=_env_int("MESON_FEEDBACK_P_NUM", 8),
            t_num=_env_int("MESON_FEEDBACK_T_NUM", 4),
            quark_residual_norm_max=_env_float("MESON_FEEDBACK_QUARK_RESIDUAL_MAX", 1e-5),
            quark_iterations=_env_int("MESON_FEEDBACK_QUARK_ITERATIONS", 200),
            gap_residual_norm_max=_env_float("MESON_FEEDBACK_GAP_RESIDUAL_MAX", 1e-5),
            outer_residual_tolerance=_env_float("MESON_FEEDBACK_OUTER_RESIDUAL_TOL", 2e-3),
            outer_finite_difference_step=_env_float("MESON_FEEDBACK_FD_STEP", 5e-3),
            outer_maximum_step=_env_float("MESON_FEEDBACK_MAXIMUM_STEP", 0.25),
            outer_max_iterations=_env_int("MESON_FEEDBACK_MAX_ITERATIONS", 12),
            outer_max_evaluations=_env_int("MESON_FEEDBACK_MAX_EVALUATIONS", 50),
        )
        push!(rows, _row(point, output, initial_source))

        if output.baseline.converged
            previous_solution = output.baseline.solution
        end
        if output.feedback !== nothing && output.feedback.result.converged
            previous_feedback = (
                mu_Q=output.feedback.result.mu_Q,
                mu_S=output.feedback.result.mu_S,
            )
        else
            previous_feedback = nothing
        end
        @printf("[freezeout-partial-feedback] sqrt(s)=%.3f GeV status=%s time=%.3fs\n",
                point.sqrt_s_NN_GeV, rows[end].status, rows[end].total_elapsed_s)
    end

    output_dir = get(ENV, "MESON_FEEDBACK_SCAN_OUTPUT_DIR", joinpath(
        PROJECT_ROOT,
        "data", "outputs", "results", "relaxtime", "analysis",
        "meson_conserved_charge_feedback_freezeout_sparse",
    ))
    csv_path = joinpath(output_dir, "partial_feedback_freezeout_scan.csv")
    summary_path = joinpath(output_dir, "README.md")
    _write_csv(csv_path, sort(rows; by=r -> r.sqrt_s_NN_GeV))
    _write_summary(summary_path, rows, median_s, median_source, grid, settings)
    println("[freezeout-partial-feedback] csv=$(csv_path)")
    println("[freezeout-partial-feedback] summary=$(summary_path)")
    return (rows=rows, csv=csv_path, summary=summary_path, grid=grid)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
