"""
Single-point partial-feedback spike for charged pion/kaon conserved charges.

For every outer `(mu_Q,mu_S)` candidate the script re-solves the five PNJL
mean fields at fixed `(T,mu_B)`, evaluates quark net B/Q/S densities, then adds
current-BU `pi+`, `pi-`, `K+`, and `K-` densities. Mesons do not enter the
mean-field stationarity equations, so every result remains diagnostic.
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(@__DIR__, "meson_conserved_charge_feedback_utils.jl"))
include(joinpath(@__DIR__, "meson_conserved_charge_feedback_runtime.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .Models
using .MesonConservedChargeFeedbackUtils
using .MesonConservedChargeFeedbackRuntime: FeedbackSettings,
                                           build_candidate_evaluator,
                                           solve_feedback_level

const DEFAULT_OUTPUT = joinpath(
    PROJECT_ROOT,
    "data",
    "outputs",
    "results",
    "relaxtime",
    "analysis",
    "meson_conserved_charge_outer_feedback_spike",
    "fixed_point_feedback.csv",
)

@inline _env_float(name::AbstractString, default::Real) = parse(Float64, get(ENV, name, string(default)))
@inline _env_int(name::AbstractString, default::Integer) = parse(Int, get(ENV, name, string(default)))
@inline _env_bool(name::AbstractString, default::Bool) = lowercase(strip(get(ENV, name, string(default)))) in ("1", "true", "yes", "on")
@inline _env_optional_float(name::AbstractString) = haskey(ENV, name) ? parse(Float64, ENV[name]) : nothing

function _settings(; refined::Bool=false)
    q_nodes = _env_int("MESON_FEEDBACK_Q_NODES", 4)
    omega_nodes = _env_int("MESON_FEEDBACK_OMEGA_NODES", 6)
    refinement = refined ? _env_int("MESON_FEEDBACK_REFINEMENT_FACTOR", 2) : 1
    return FeedbackSettings(
        label=refined ? "refined" : "coarse",
        qmax=_env_float("MESON_FEEDBACK_QMAX", 4.0),
        q_nodes=q_nodes * refinement,
        omega_min=_env_float("MESON_FEEDBACK_OMEGA_MIN", 0.05),
        omega_max=_env_float("MESON_FEEDBACK_OMEGA_MAX", 3.0),
        omega_nodes=omega_nodes * refinement,
        eta=_env_float("MESON_FEEDBACK_ETA", 1e-6),
        density_policy=:x_min_cut,
        bose_x_min=_env_float("MESON_FEEDBACK_BOSE_X_MIN", 0.05),
    )
end

function _candidate_evaluator(
    model,
    T_fm::Float64,
    muB_fm::Float64,
    branch_seed::Vector{Float64},
    settings;
    p_num::Int,
    t_num::Int,
    gap_residual_norm_max::Float64,
)
    return build_candidate_evaluator(
        model,
        T_fm,
        muB_fm,
        branch_seed,
        settings;
        p_num=p_num,
        t_num=t_num,
        gap_residual_norm_max=gap_residual_norm_max,
    )
end

function _solve_level(
    model,
    baseline,
    T_fm::Float64,
    muB_fm::Float64,
    initial_mu_Q::Float64,
    initial_mu_S::Float64,
    settings;
    baseline_mu_Q::Union{Nothing,Float64}=nothing,
    baseline_mu_S::Union{Nothing,Float64}=nothing,
    p_num::Int,
    t_num::Int,
    rho0::Float64,
    target_ratio::Float64,
    rho_S_target::Float64,
)
    return solve_feedback_level(
        model,
        baseline,
        T_fm,
        muB_fm,
        initial_mu_Q,
        initial_mu_S,
        settings;
        baseline_mu_Q=baseline_mu_Q,
        baseline_mu_S=baseline_mu_S,
        p_num=p_num,
        t_num=t_num,
        charge_to_baryon_ratio=target_ratio,
        rho0=rho0,
        rho_S_target=rho_S_target,
        gap_residual_norm_max=_env_float("MESON_FEEDBACK_GAP_RESIDUAL_MAX", 1e-5),
        residual_tolerance=_env_float("MESON_FEEDBACK_OUTER_RESIDUAL_TOL", 2e-3),
        finite_difference_step=_env_float("MESON_FEEDBACK_FD_STEP", 5e-3),
        maximum_step=_env_float("MESON_FEEDBACK_MAXIMUM_STEP", 0.25),
        max_iterations=_env_int("MESON_FEEDBACK_MAX_ITERATIONS", 6),
        max_evaluations=_env_int("MESON_FEEDBACK_MAX_EVALUATIONS", 30),
    )
end

@inline _ratio(numerator::Real, denominator::Real) = denominator > 0.0 ? numerator / denominator : NaN

function _row(stage::String, level, payload, mu_Q::Float64, mu_S::Float64, residual;
              converged::Bool, reason::Symbol, iterations::Int, evaluations::Int, verdict::String="")
    meson = charged_meson_conserved_densities(payload)
    total = total_conserved_charge_residual(
        (rho_B=payload.rho_B_q, rho_Q=payload.rho_Q_q, rho_S=payload.rho_S_q),
        meson;
        charge_to_baryon_ratio=_env_float("MESON_FEEDBACK_TARGET_QB", 0.4),
        strangeness_density_target=_env_float("MESON_FEEDBACK_RHOS_TARGET", 0.0),
        rho0=_env_float("MESON_FEEDBACK_RHO0", 0.16),
    )
    return (
        stage=stage,
        numeric_level=level.settings.label,
        converged=converged,
        reason=String(reason),
        verdict=verdict,
        T_MeV=_env_float("MESON_FEEDBACK_T_MEV", 170.0),
        muB_MeV=_env_float("MESON_FEEDBACK_MUB_MEV", 240.0),
        muQ_MeV=mu_Q * ħc_MeV_fm,
        muS_MeV=mu_S * ħc_MeV_fm,
        mu_u_MeV=payload.mu_u * ħc_MeV_fm,
        mu_d_MeV=payload.mu_d * ħc_MeV_fm,
        mu_s_MeV=payload.mu_s * ħc_MeV_fm,
        rho_B_q=payload.rho_B_q,
        rho_Q_q=payload.rho_Q_q,
        rho_S_q=payload.rho_S_q,
        rho_Q_M=meson.rho_Q,
        rho_S_M=meson.rho_S,
        rho_Q_total=total.rho_Q_total,
        rho_S_total=total.rho_S_total,
        charge_residual=total.normalized[1],
        strangeness_residual=total.normalized[2],
        residual_norm=total.norm,
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
        gap_residual_norm=payload.gap_residual_norm,
        outer_iterations=iterations,
        outer_evaluations=evaluations,
        outer_initial_muQ_MeV=level.initial_mu_Q * ħc_MeV_fm,
        outer_initial_muS_MeV=level.initial_mu_S * ħc_MeV_fm,
        unique_candidate_count=level.unique_candidate_count,
        scheme="current",
        density_policy=String(level.settings.density_policy),
        bose_x_min=level.settings.bose_x_min,
        qmax=level.settings.qmax,
        q_nodes=level.settings.q_nodes,
        omega_min=level.settings.omega_min,
        omega_max=level.settings.omega_max,
        omega_nodes=level.settings.omega_nodes,
        eta=level.settings.eta,
    )
end

@inline function _relative_change(new::Real, old::Real)
    return abs(Float64(new) - Float64(old)) / max(abs(Float64(old)), eps(Float64))
end

function _verdict(coarse, refined, T_fm::Float64, initial_mu_Q::Float64, initial_mu_S::Float64)
    coarse.result.converged || return "blocked:coarse_outer_not_converged"
    refined.result.converged || return "blocked:refined_outer_not_converged"

    coarse_payload = coarse.result.payload
    refined_payload = refined.result.payload
    feedback_mu = max(
        abs(refined.result.mu_Q - initial_mu_Q),
        abs(refined.result.mu_S - initial_mu_S),
    )
    node_mu = max(
        abs(refined.result.mu_Q - coarse.result.mu_Q),
        abs(refined.result.mu_S - coarse.result.mu_S),
    )
    observables = (
        :n_pi_plus,
        :n_pi_minus,
        :n_K_plus,
        :n_K_minus,
    )
    component_stable = all(observables) do name
        baseline_coarse = getproperty(coarse.baseline_payload, name)
        baseline_refined = getproperty(refined.baseline_payload, name)
        corrected_coarse = getproperty(coarse_payload, name)
        corrected_refined = getproperty(refined_payload, name)
        feedback_change = _relative_change(corrected_refined, baseline_refined)
        node_change = max(
            _relative_change(baseline_refined, baseline_coarse),
            _relative_change(corrected_refined, corrected_coarse),
        )
        node_change <= max(0.5 * feedback_change, 0.05)
    end

    baseline_plus_coarse = _ratio(coarse.baseline_payload.n_K_plus, coarse.baseline_payload.n_pi_plus)
    baseline_minus_coarse = _ratio(coarse.baseline_payload.n_K_minus, coarse.baseline_payload.n_pi_minus)
    baseline_plus_refined = _ratio(refined.baseline_payload.n_K_plus, refined.baseline_payload.n_pi_plus)
    baseline_minus_refined = _ratio(refined.baseline_payload.n_K_minus, refined.baseline_payload.n_pi_minus)
    corrected_plus_coarse = _ratio(coarse_payload.n_K_plus, coarse_payload.n_pi_plus)
    corrected_minus_coarse = _ratio(coarse_payload.n_K_minus, coarse_payload.n_pi_minus)
    corrected_plus_refined = _ratio(refined_payload.n_K_plus, refined_payload.n_pi_plus)
    corrected_minus_refined = _ratio(refined_payload.n_K_minus, refined_payload.n_pi_minus)
    plus_effect = _relative_change(corrected_plus_refined, baseline_plus_refined)
    minus_effect = _relative_change(corrected_minus_refined, baseline_minus_refined)
    plus_node = max(
        _relative_change(baseline_plus_refined, baseline_plus_coarse),
        _relative_change(corrected_plus_refined, corrected_plus_coarse),
    )
    minus_node = max(
        _relative_change(baseline_minus_refined, baseline_minus_coarse),
        _relative_change(corrected_minus_refined, corrected_minus_coarse),
    )
    ratio_stable = plus_node <= max(0.5 * plus_effect, 0.05) &&
                   minus_node <= max(0.5 * minus_effect, 0.05)
    ratio_effect = max(plus_effect, minus_effect)
    stable = node_mu <= max(0.5 * feedback_mu, 0.01 * T_fm) &&
             component_stable && ratio_stable
    stable || return "blocked:node_refinement_not_subleading"
    significant = feedback_mu > 0.05 * T_fm || ratio_effect > 0.10
    return significant ? "full-feedback-candidate" : "skip-full-feedback"
end

function _write_csv(path::String, rows)
    isempty(rows) && return
    mkpath(dirname(path))
    columns = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(String.(columns), ','))
        for row in rows
            values = (replace(string(getproperty(row, col)), ',' => ';') for col in columns)
            println(io, join(values, ','))
        end
    end
end

function main()
    T_fm = _env_float("MESON_FEEDBACK_T_MEV", 170.0) / ħc_MeV_fm
    muB_fm = _env_float("MESON_FEEDBACK_MUB_MEV", 240.0) / ħc_MeV_fm
    target_ratio = _env_float("MESON_FEEDBACK_TARGET_QB", 0.4)
    rho_S_target = _env_float("MESON_FEEDBACK_RHOS_TARGET", 0.0)
    rho0 = _env_float("MESON_FEEDBACK_RHO0", 0.16)
    p_num = _env_int("MESON_FEEDBACK_P_NUM", 8)
    t_num = _env_int("MESON_FEEDBACK_T_NUM", 4)
    model = create_model(:PNJL)
    mode = FixedMuBConservedCharges(muB_fm, target_ratio, rho_S_target)
    baseline = solve(
        model,
        mode,
        T_fm;
        p_num=p_num,
        t_num=t_num,
        residual_norm_max=_env_float("MESON_FEEDBACK_QUARK_RESIDUAL_MAX", 1e-5),
        iterations=_env_int("MESON_FEEDBACK_QUARK_ITERATIONS", 200),
    )
    baseline.converged || error("quark-only baseline did not converge: $(baseline.residual_norm)")
    conserved_mu = conserved_mu_from_flavor(baseline.mu_vec...)
    initial_muQ_MeV = _env_optional_float("MESON_FEEDBACK_INITIAL_MUQ_MEV")
    initial_muS_MeV = _env_optional_float("MESON_FEEDBACK_INITIAL_MUS_MEV")
    outer_initial_mu_Q = initial_muQ_MeV === nothing ? conserved_mu.mu_Q : initial_muQ_MeV / ħc_MeV_fm
    outer_initial_mu_S = initial_muS_MeV === nothing ? conserved_mu.mu_S : initial_muS_MeV / ħc_MeV_fm

    coarse = _solve_level(
        model,
        baseline,
        T_fm,
        muB_fm,
        outer_initial_mu_Q,
        outer_initial_mu_S,
        _settings();
        baseline_mu_Q=conserved_mu.mu_Q,
        baseline_mu_S=conserved_mu.mu_S,
        p_num=p_num,
        t_num=t_num,
        rho0=rho0,
        target_ratio=target_ratio,
        rho_S_target=rho_S_target,
    )
    run_refined = _env_bool("MESON_FEEDBACK_RUN_REFINED_OUTER", true)
    refined = run_refined ? _solve_level(
        model,
        baseline,
        T_fm,
        muB_fm,
        coarse.result.converged ? coarse.result.mu_Q : conserved_mu.mu_Q,
        coarse.result.converged ? coarse.result.mu_S : conserved_mu.mu_S,
        _settings(; refined=true);
        baseline_mu_Q=conserved_mu.mu_Q,
        baseline_mu_S=conserved_mu.mu_S,
        p_num=p_num,
        t_num=t_num,
        rho0=rho0,
        target_ratio=target_ratio,
        rho_S_target=rho_S_target,
    ) : nothing
    verdict = refined === nothing ? "blocked:refined_outer_not_run" : _verdict(
        coarse,
        refined,
        T_fm,
        conserved_mu.mu_Q,
        conserved_mu.mu_S,
    )

    rows = NamedTuple[]
    coarse_baseline_residual = total_conserved_charge_residual(
        (rho_B=coarse.baseline_payload.rho_B_q, rho_Q=coarse.baseline_payload.rho_Q_q, rho_S=coarse.baseline_payload.rho_S_q),
        charged_meson_conserved_densities(coarse.baseline_payload);
        charge_to_baryon_ratio=target_ratio,
        strangeness_density_target=rho_S_target,
        rho0=rho0,
    )
    push!(rows, _row(
        "quark_only_seed",
        coarse,
        coarse.baseline_payload,
        conserved_mu.mu_Q,
        conserved_mu.mu_S,
        coarse_baseline_residual;
        converged=baseline.converged,
        reason=:quark_only_seed,
        iterations=0,
        evaluations=1,
        verdict=verdict,
    ))
    push!(rows, _row(
        "partial_feedback",
        coarse,
        coarse.result.payload,
        coarse.result.mu_Q,
        coarse.result.mu_S,
        coarse.result.residual;
        converged=coarse.result.converged,
        reason=coarse.result.reason,
        iterations=coarse.result.iterations,
        evaluations=coarse.result.evaluation_count,
        verdict=verdict,
    ))

    if refined !== nothing
        push!(rows, _row(
            "quark_only_seed",
            refined,
            refined.baseline_payload,
            conserved_mu.mu_Q,
            conserved_mu.mu_S,
            nothing;
            converged=baseline.converged,
            reason=:refined_seed,
            iterations=0,
            evaluations=1,
            verdict=verdict,
        ))
        push!(rows, _row(
            "partial_feedback",
            refined,
            refined.result.payload,
            refined.result.mu_Q,
            refined.result.mu_S,
            refined.result.residual;
            converged=refined.result.converged,
            reason=refined.result.reason,
            iterations=refined.result.iterations,
            evaluations=refined.result.evaluation_count,
            verdict=verdict,
        ))
    end

    output = get(ENV, "MESON_FEEDBACK_OUTPUT", DEFAULT_OUTPUT)
    _write_csv(output, rows)
    println("[meson-feedback] output=$(output)")
    println("[meson-feedback] verdict=$(verdict)")
    println("[meson-feedback] coarse_converged=$(coarse.result.converged) residual=$(coarse.result.residual_norm) evaluations=$(coarse.result.evaluation_count)")
    refined === nothing || println("[meson-feedback] refined_converged=$(refined.result.converged) residual=$(refined.result.residual_norm) evaluations=$(refined.result.evaluation_count)")
    return (rows=rows, coarse=coarse, refined=refined, verdict=verdict, output=output)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
