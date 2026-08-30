"""
Hot-start benchmark for one charged-meson partial-feedback point.

The first complete point is a warm-up and is excluded from the reported
statistics.  Every measured sample creates a fresh evaluator/cache so the
outer Newton loop cannot obtain zero-cost hits from a previous sample.
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "analysis", "relaxtime", "meson_conserved_charge_feedback_utils.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "analysis", "relaxtime", "meson_conserved_charge_feedback_runtime.jl"))

using JSON
using Statistics: mean, median
using Printf: @printf
using .Constants_PNJL: ħc_MeV_fm
using .Models
using .MesonConservedChargeFeedbackUtils
using .MesonConservedChargeFeedbackRuntime: FeedbackSettings, solve_partial_feedback_point

@inline _env_float(name::AbstractString, default::Real) = parse(Float64, get(ENV, name, string(default)))
@inline _env_int(name::AbstractString, default::Integer) = parse(Int, get(ENV, name, string(default)))

function _settings()
    return FeedbackSettings(
        label="benchmark_coarse",
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

function _run_point(settings, T_fm::Float64, muB_fm::Float64)
    model = create_model(:PNJL)
    return solve_partial_feedback_point(
        model,
        T_fm,
        muB_fm,
        settings;
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
end

@inline _ratio(numerator::Real, denominator::Real) = denominator > 0.0 ? numerator / denominator : NaN

function _sample_row(sample::Int, elapsed_s::Float64, output)
    feedback = output.feedback
    if feedback === nothing
        return (
            sample=sample,
            elapsed_s=elapsed_s,
            baseline_elapsed_s=output.baseline_elapsed_s,
            outer_elapsed_s=NaN,
            candidate_elapsed_s=NaN,
            gap_elapsed_s=NaN,
            density_elapsed_s=NaN,
            unique_candidate_count=0,
            outer_evaluations=0,
            outer_iterations=0,
            converged=false,
            residual_norm=Inf,
            charge_residual=NaN,
            strangeness_residual=NaN,
            Kplus_over_piplus=NaN,
            Kminus_over_piminus=NaN,
            reason=String(output.reason),
        )
    end
    result = feedback.result
    payload = result.payload
    timing = feedback.timing
    residual = result.residual
    return (
        sample=sample,
        elapsed_s=elapsed_s,
        baseline_elapsed_s=output.baseline_elapsed_s,
        outer_elapsed_s=feedback.outer_elapsed_s,
        candidate_elapsed_s=timing.candidate_elapsed_s,
        gap_elapsed_s=timing.gap_elapsed_s,
        density_elapsed_s=timing.density_elapsed_s,
        unique_candidate_count=timing.unique_candidate_count,
        outer_evaluations=result.evaluation_count,
        outer_iterations=result.iterations,
        converged=result.converged,
        residual_norm=result.residual_norm,
        charge_residual=residual.normalized[1],
        strangeness_residual=residual.normalized[2],
        Kplus_over_piplus=_ratio(payload.n_K_plus, payload.n_pi_plus),
        Kminus_over_piminus=_ratio(payload.n_K_minus, payload.n_pi_minus),
        reason=String(result.reason),
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

function _summary_dict(rows, warmup_elapsed_s::Float64, settings, T_MeV::Float64, muB_MeV::Float64)
    times = Float64[row.elapsed_s for row in rows]
    med = median(times)
    return Dict(
        "kind" => "diagnostic_hot_start_partial_feedback",
        "T_MeV" => T_MeV,
        "muB_MeV" => muB_MeV,
        "warmup_elapsed_s" => warmup_elapsed_s,
        "sample_count" => length(rows),
        "elapsed_s" => Dict(
            "min" => minimum(times),
            "median" => med,
            "mean" => mean(times),
            "max" => maximum(times),
        ),
        "recommended_sqrts_grid_GeV" => choose_freezeout_sqrts_grid(med),
        "settings" => Dict(
            "p_num" => _env_int("MESON_FEEDBACK_P_NUM", 8),
            "t_num" => _env_int("MESON_FEEDBACK_T_NUM", 4),
            "qmax" => settings.qmax,
            "q_nodes" => settings.q_nodes,
            "omega_min" => settings.omega_min,
            "omega_max" => settings.omega_max,
            "omega_nodes" => settings.omega_nodes,
            "density_policy" => String(settings.density_policy),
            "bose_x_min" => settings.bose_x_min,
            "outer_max_iterations" => _env_int("MESON_FEEDBACK_MAX_ITERATIONS", 12),
            "outer_max_evaluations" => _env_int("MESON_FEEDBACK_MAX_EVALUATIONS", 50),
        ),
        "samples" => [Dict(String(k) => getproperty(row, k) for k in keys(row)) for row in rows],
    )
end

function main()
    T_MeV = _env_float("MESON_FEEDBACK_BENCH_T_MEV", 170.0)
    muB_MeV = _env_float("MESON_FEEDBACK_BENCH_MUB_MEV", 240.0)
    samples = _env_int("MESON_FEEDBACK_BENCH_SAMPLES", 3)
    samples > 0 || throw(ArgumentError("MESON_FEEDBACK_BENCH_SAMPLES must be positive"))
    settings = _settings()
    T_fm = T_MeV / ħc_MeV_fm
    muB_fm = muB_MeV / ħc_MeV_fm

    GC.gc()
    warmup_start = time_ns()
    warmup = _run_point(settings, T_fm, muB_fm)
    warmup_elapsed_s = (time_ns() - warmup_start) / 1.0e9
    if !warmup.converged
        residual_text = warmup.feedback === nothing ? "NaN" : string(warmup.feedback.result.residual_norm)
        history_text = warmup.feedback === nothing ? "[]" : string(warmup.feedback.result.history)
        @printf("[partial-feedback-benchmark] warmup status=%s residual=%s history=%s\n", warmup.reason, residual_text, history_text)
    end

    rows = NamedTuple[]
    for sample in 1:samples
        GC.gc()
        t0 = time_ns()
        output = _run_point(settings, T_fm, muB_fm)
        elapsed_s = (time_ns() - t0) / 1.0e9
        push!(rows, _sample_row(sample, elapsed_s, output))
    end

    output_dir = get(ENV, "MESON_FEEDBACK_BENCH_OUTPUT_DIR", joinpath(
        PROJECT_ROOT,
        "data", "outputs", "perf", "relaxtime", "meson_conserved_charge_partial_feedback",
    ))
    csv_path = joinpath(output_dir, "hot_start_samples.csv")
    json_path = joinpath(output_dir, "hot_start_summary.json")
    _write_csv(csv_path, rows)
    summary = _summary_dict(rows, warmup_elapsed_s, settings, T_MeV, muB_MeV)
    mkpath(dirname(json_path))
    open(json_path, "w") do io
        print(io, JSON.json(summary; pretty=true))
    end

    med = summary["elapsed_s"]["median"]
    recommended = summary["recommended_sqrts_grid_GeV"]
    @printf("[partial-feedback-benchmark] warmup=%.3fs samples=%d median=%.3fs mean=%.3fs\n",
            warmup_elapsed_s, samples, med, summary["elapsed_s"]["mean"])
    println("[partial-feedback-benchmark] recommended_sqrts=$(recommended)")
    println("[partial-feedback-benchmark] csv=$(csv_path)")
    println("[partial-feedback-benchmark] summary=$(json_path)")
    return (rows=rows, summary=summary, csv=csv_path, json=json_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
