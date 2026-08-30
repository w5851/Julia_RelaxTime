"""
Internal runtime shared by the charged-meson conserved-charge diagnostics.

This is intentionally an analysis-layer module.  It reuses the stable Models
solver and the current BU density kernel, but it is not a public workflow API
and it must not be used to label a result as a thermodynamically self-consistent
meson-feedback equilibrium.
"""
module MesonConservedChargeFeedbackRuntime

using LinearAlgebra: norm

using Main.Constants_PNJL: G_fm2, K_fm5, ħc_MeV_fm
using Main.Models
using Main.MesonConservedChargeFeedbackUtils
using Main.MesonInteractionKernel: build_full_kmt_interaction, charged_coupling
using Main.MesonDensity: phase_shift_meson_number_density

export FeedbackSettings
export build_candidate_evaluator
export solve_feedback_level
export solve_partial_feedback_point
export solve_quark_only_bu_ab_point
export candidate_timing_summary

Base.@kwdef struct FeedbackSettings
    label::String = "coarse"
    qmax::Float64 = 4.0
    q_nodes::Int = 4
    omega_min::Float64 = 0.05
    omega_max::Float64 = 3.0
    omega_nodes::Int = 8
    eta::Float64 = 1e-6
    density_policy::Symbol = :x_min_cut
    bose_x_min::Float64 = 0.05
end

@inline function _finite(name::AbstractString, value::Real)::Float64
    x = Float64(value)
    isfinite(x) || throw(ArgumentError("$(name) must be finite, got $(value)"))
    return x
end

function _validate_settings(settings::FeedbackSettings)
    settings.qmax > 0.0 || throw(ArgumentError("qmax must be positive"))
    settings.q_nodes > 0 || throw(ArgumentError("q_nodes must be positive"))
    settings.omega_min > 0.0 || throw(ArgumentError("omega_min must be positive"))
    settings.omega_max > settings.omega_min || throw(ArgumentError("omega_max must exceed omega_min"))
    settings.omega_nodes > 0 || throw(ArgumentError("omega_nodes must be positive"))
    settings.eta > 0.0 || throw(ArgumentError("eta must be positive"))
    settings.bose_x_min > 0.0 || throw(ArgumentError("bose_x_min must be positive"))
    return nothing
end

@inline function _phase_density(
    meson::Symbol,
    mu_M::Float64,
    qp,
    tp,
    settings::FeedbackSettings;
    interaction=nothing,
)
    return phase_shift_meson_number_density(
        meson,
        qp,
        tp;
        degeneracy=1,
        μ=mu_M,
        scheme=:current,
        qmax=settings.qmax,
        q_nodes=settings.q_nodes,
        omega_min=settings.omega_min,
        omega_max=settings.omega_max,
        omega_nodes=settings.omega_nodes,
        eta=settings.eta,
        density_policy=settings.density_policy,
        bose_x_min=settings.bose_x_min,
        interaction=interaction,
    )
end

function _model_kmt_couplings(model)
    if hasproperty(model, :params)
        params = getproperty(model, :params)
        if hasproperty(params, :G_fm2) && hasproperty(params, :K_fm5)
            return (G=Float64(params.G_fm2), K=Float64(params.K_fm5), source=:model_params)
        end
    end
    return (G=Float64(G_fm2), K=Float64(K_fm5), source=:constants_fallback)
end

function _charged_density_set(meson_mu, qp, tp, settings::FeedbackSettings; interaction=nothing)
    density_start = time_ns()
    pi_plus = _phase_density(:pi_plus, meson_mu.mu_pi_plus, qp, tp, settings; interaction=interaction)
    pi_minus = _phase_density(:pi_minus, meson_mu.mu_pi_minus, qp, tp, settings; interaction=interaction)
    K_plus = _phase_density(:K_plus, meson_mu.mu_K_plus, qp, tp, settings; interaction=interaction)
    K_minus = _phase_density(:K_minus, meson_mu.mu_K_minus, qp, tp, settings; interaction=interaction)
    return (
        pi_plus=pi_plus,
        pi_minus=pi_minus,
        K_plus=K_plus,
        K_minus=K_minus,
        elapsed_s=(time_ns() - density_start) / 1.0e9,
    )
end

function build_candidate_evaluator(
    model,
    T_fm::Float64,
    muB_fm::Float64,
    branch_seed::Vector{Float64},
    settings::FeedbackSettings;
    p_num::Int,
    t_num::Int,
    gap_residual_norm_max::Float64,
)
    _validate_settings(settings)
    cache = Dict{Tuple{Float64,Float64},NamedTuple}()

    function evaluate(mu_Q::Real, mu_S::Real)
        key = (round(Float64(mu_Q); digits=12), round(Float64(mu_S); digits=12))
        return get!(cache, key) do
            candidate_start = time_ns()
            flavor = flavor_mu_from_bqs(muB_fm, key[1], key[2])
            mu_vec = Float64[flavor.mu_u, flavor.mu_d, flavor.mu_s]

            gap_start = time_ns()
            state = solve_gap(
                model,
                T_fm,
                mu_vec;
                solver_backend=:models,
                initial_guess=branch_seed,
                residual_norm_max=gap_residual_norm_max,
                xi=0.0,
                p_num=p_num,
                t_num=t_num,
            )
            gap_elapsed_s = (time_ns() - gap_start) / 1.0e9
            x_state = collect(state_vector(state))
            gap_norm = norm(gap_residual(
                model,
                x_state,
                T_fm,
                mu_vec;
                xi=0.0,
                p_num=p_num,
                t_num=t_num,
            ))
            isfinite(gap_norm) && gap_norm <= gap_residual_norm_max || throw(ArgumentError(
                "candidate gap residual $(gap_norm) exceeds $(gap_residual_norm_max)",
            ))

            masses = calculate_mass_vec(model, meanfield_state(x_state).phi)
            qp = (
                m=(u=Float64(masses[1]), d=Float64(masses[2]), s=Float64(masses[3])),
                μ=(u=mu_vec[1], d=mu_vec[2], s=mu_vec[3]),
                # Preserve the solved condensates for the phi-native legacy
                # coupling adapter; A is still built for the polarization.
                phi=(u=Float64(x_state[1]), d=Float64(x_state[2]), s=Float64(x_state[3])),
            )
            tp = (
                T=T_fm,
                Φ=Float64(x_state[4]),
                Φbar=Float64(x_state[5]),
                ξ=0.0,
            )

            quark_rho = model_rho(
                model,
                x_state,
                mu_vec,
                T_fm;
                xi=0.0,
                p_num=p_num,
                t_num=t_num,
            )
            quark_bqs = conserved_densities_from_flavor(quark_rho)
            meson_mu = charged_meson_chemical_potentials(key[1], key[2])

            density_start = time_ns()
            pi_plus = _phase_density(:pi_plus, meson_mu.mu_pi_plus, qp, tp, settings)
            pi_minus = _phase_density(:pi_minus, meson_mu.mu_pi_minus, qp, tp, settings)
            K_plus = _phase_density(:K_plus, meson_mu.mu_K_plus, qp, tp, settings)
            K_minus = _phase_density(:K_minus, meson_mu.mu_K_minus, qp, tp, settings)
            density_elapsed_s = (time_ns() - density_start) / 1.0e9

            return (
                mu_u=mu_vec[1],
                mu_d=mu_vec[2],
                mu_s=mu_vec[3],
                x_state=x_state,
                gap_residual_norm=gap_norm,
                rho_B_q=quark_bqs.rho_B,
                rho_Q_q=quark_bqs.rho_Q,
                rho_S_q=quark_bqs.rho_S,
                n_pi_plus=Float64(pi_plus.density),
                n_pi_minus=Float64(pi_minus.density),
                n_K_plus=Float64(K_plus.density),
                n_K_minus=Float64(K_minus.density),
                pi_plus_status=pi_plus.status,
                pi_minus_status=pi_minus.status,
                K_plus_status=K_plus.status,
                K_minus_status=K_minus.status,
                pi_plus_min_E_minus_mu=pi_plus.min_E_minus_mu,
                pi_minus_min_E_minus_mu=pi_minus.min_E_minus_mu,
                K_plus_min_E_minus_mu=K_plus.min_E_minus_mu,
                K_minus_min_E_minus_mu=K_minus.min_E_minus_mu,
                gap_elapsed_s=gap_elapsed_s,
                density_elapsed_s=density_elapsed_s,
                candidate_elapsed_s=(time_ns() - candidate_start) / 1.0e9,
            )
        end
    end

    return evaluate, cache
end

function candidate_timing_summary(cache)
    isempty(cache) && return (
        unique_candidate_count=0,
        gap_elapsed_s=0.0,
        density_elapsed_s=0.0,
        candidate_elapsed_s=0.0,
        gap_mean_ms=NaN,
        density_mean_ms=NaN,
    )
    gap_s = sum(Float64(getproperty(v, :gap_elapsed_s)) for v in values(cache))
    density_s = sum(Float64(getproperty(v, :density_elapsed_s)) for v in values(cache))
    n = length(cache)
    return (
        unique_candidate_count=n,
        gap_elapsed_s=gap_s,
        density_elapsed_s=density_s,
    candidate_elapsed_s=sum(Float64(getproperty(v, :candidate_elapsed_s)) for v in values(cache)),
        gap_mean_ms=1000.0 * gap_s / n,
        density_mean_ms=1000.0 * density_s / n,
    )
end

function solve_feedback_level(
    model,
    baseline,
    T_fm::Float64,
    muB_fm::Float64,
    initial_mu_Q::Float64,
    initial_mu_S::Float64,
    settings::FeedbackSettings;
    baseline_mu_Q::Union{Nothing,Float64}=nothing,
    baseline_mu_S::Union{Nothing,Float64}=nothing,
    p_num::Int,
    t_num::Int,
    rho0::Float64,
    target_ratio::Float64,
    rho_S_target::Float64,
    gap_residual_norm_max::Float64=1e-5,
    residual_tolerance::Float64=2e-3,
    finite_difference_step::Float64=5e-3,
    maximum_step::Float64=0.25,
    max_iterations::Int=6,
    max_evaluations::Int=30,
)
    evaluator, cache = build_candidate_evaluator(
        model,
        T_fm,
        muB_fm,
        collect(baseline.x_state),
        settings;
        p_num=p_num,
        t_num=t_num,
        gap_residual_norm_max=gap_residual_norm_max,
    )
    baseline_Q = baseline_mu_Q === nothing ? initial_mu_Q : baseline_mu_Q
    baseline_S = baseline_mu_S === nothing ? initial_mu_S : baseline_mu_S
    baseline_payload = evaluator(baseline_Q, baseline_S)

    outer_start = time_ns()
    result = solve_outer_conserved_charge_feedback(
        evaluator,
        initial_mu_Q,
        initial_mu_S;
        charge_to_baryon_ratio=target_ratio,
        strangeness_density_target=rho_S_target,
        rho0=rho0,
        residual_tolerance=residual_tolerance,
        finite_difference_step=finite_difference_step,
        maximum_step=maximum_step,
        max_iterations=max_iterations,
        max_evaluations=max_evaluations,
    )
    outer_elapsed_s = (time_ns() - outer_start) / 1.0e9
    timing = candidate_timing_summary(cache)
    return (
        settings=settings,
        initial_mu_Q=initial_mu_Q,
        initial_mu_S=initial_mu_S,
        baseline_payload=baseline_payload,
        result=result,
        cache=cache,
        outer_elapsed_s=outer_elapsed_s,
        timing=timing,
        unique_candidate_count=timing.unique_candidate_count,
    )
end

function solve_partial_feedback_point(
    model,
    T_fm::Float64,
    muB_fm::Float64,
    settings::FeedbackSettings;
    baseline_seed=nothing,
    initial_mu_Q::Union{Nothing,Float64}=nothing,
    initial_mu_S::Union{Nothing,Float64}=nothing,
    target_ratio::Float64=0.4,
    rho_S_target::Float64=0.0,
    rho0::Float64=0.16,
    p_num::Int=8,
    t_num::Int=4,
    quark_residual_norm_max::Float64=1e-5,
    quark_iterations::Int=200,
    gap_residual_norm_max::Float64=1e-5,
    outer_residual_tolerance::Float64=2e-3,
    outer_finite_difference_step::Float64=5e-3,
    outer_maximum_step::Float64=0.25,
    outer_max_iterations::Int=6,
    outer_max_evaluations::Int=30,
)
    mode = FixedMuBConservedCharges(muB_fm, target_ratio, rho_S_target)
    baseline_kwargs = (
        p_num=p_num,
        t_num=t_num,
        residual_norm_max=quark_residual_norm_max,
        iterations=quark_iterations,
        xi=0.0,
    )
    baseline_kwargs = baseline_seed === nothing ? baseline_kwargs : (; baseline_kwargs..., seed_guess=baseline_seed)

    baseline_start = time_ns()
    baseline = solve(model, mode, T_fm; baseline_kwargs...)
    baseline_elapsed_s = (time_ns() - baseline_start) / 1.0e9
    baseline.converged || return (
        converged=false,
        reason=:quark_only_not_converged,
        message="quark-only baseline residual $(baseline.residual_norm) exceeds $(quark_residual_norm_max)",
        baseline=baseline,
        baseline_elapsed_s=baseline_elapsed_s,
        feedback=nothing,
        total_elapsed_s=baseline_elapsed_s,
    )

    conserved_mu = conserved_mu_from_flavor(baseline.mu_vec...)
    baseline_rho = model_rho(
        model,
        baseline.x_state,
        baseline.mu_vec,
        T_fm;
        xi=0.0,
        p_num=p_num,
        t_num=t_num,
    )
    baseline_bqs = conserved_densities_from_flavor(baseline_rho)
    start_Q = initial_mu_Q === nothing ? conserved_mu.mu_Q : Float64(initial_mu_Q)
    start_S = initial_mu_S === nothing ? conserved_mu.mu_S : Float64(initial_mu_S)
    feedback = solve_feedback_level(
        model,
        baseline,
        T_fm,
        muB_fm,
        start_Q,
        start_S,
        settings;
        baseline_mu_Q=conserved_mu.mu_Q,
        baseline_mu_S=conserved_mu.mu_S,
        p_num=p_num,
        t_num=t_num,
        rho0=rho0,
        target_ratio=target_ratio,
        rho_S_target=rho_S_target,
        gap_residual_norm_max=gap_residual_norm_max,
        residual_tolerance=outer_residual_tolerance,
        finite_difference_step=outer_finite_difference_step,
        maximum_step=outer_maximum_step,
        max_iterations=outer_max_iterations,
        max_evaluations=outer_max_evaluations,
    )
    return (
        converged=feedback.result.converged,
        reason=feedback.result.reason,
        message=feedback.result.message,
        baseline=baseline,
        baseline_mu=conserved_mu,
        baseline_bqs=baseline_bqs,
        baseline_elapsed_s=baseline_elapsed_s,
        feedback=feedback,
        total_elapsed_s=baseline_elapsed_s + feedback.outer_elapsed_s,
    )
end

"""
    solve_quark_only_bu_ab_point(model, T_fm, muB_fm, settings; kwargs...)

Solve the existing BQS quark-only equilibrium once and evaluate the four
charged BU densities twice on that same state: first with the legacy effective
couplings and then with the full charged KMT couplings `K12`/`K45`.  The
calculation is intentionally post-processing only.  In particular, the full
kernel does not alter `Omega` or the BQS stationarity equations here.
"""
function solve_quark_only_bu_ab_point(
    model,
    T_fm::Float64,
    muB_fm::Float64,
    settings::FeedbackSettings;
    baseline_seed=nothing,
    target_ratio::Float64=0.4,
    rho_S_target::Float64=0.0,
    rho0::Float64=0.16,
    p_num::Int=8,
    t_num::Int=4,
    quark_residual_norm_max::Float64=1e-5,
    quark_iterations::Int=200,
)
    _validate_settings(settings)
    mode = FixedMuBConservedCharges(muB_fm, target_ratio, rho_S_target)
    baseline_kwargs = (
        p_num=p_num,
        t_num=t_num,
        residual_norm_max=quark_residual_norm_max,
        iterations=quark_iterations,
        xi=0.0,
    )
    baseline_kwargs = baseline_seed === nothing ? baseline_kwargs : (; baseline_kwargs..., seed_guess=baseline_seed)

    baseline_start = time_ns()
    baseline = solve(model, mode, T_fm; baseline_kwargs...)
    baseline_elapsed_s = (time_ns() - baseline_start) / 1.0e9
    if !baseline.converged
        return (
            converged=false,
            reason=:quark_only_not_converged,
            message="quark-only baseline residual $(baseline.residual_norm) exceeds $(quark_residual_norm_max)",
            baseline=baseline,
            baseline_mu=nothing,
            baseline_bqs=nothing,
            meson_mu=nothing,
            legacy=nothing,
            full=nothing,
            full_kernel=nothing,
            kernel_couplings=nothing,
            baseline_elapsed_s=baseline_elapsed_s,
            legacy_density_elapsed_s=NaN,
            full_density_elapsed_s=NaN,
            total_elapsed_s=baseline_elapsed_s,
        )
    end

    conserved_mu = conserved_mu_from_flavor(baseline.mu_vec...)
    baseline_rho = model_rho(
        model,
        baseline.x_state,
        baseline.mu_vec,
        T_fm;
        xi=0.0,
        p_num=p_num,
        t_num=t_num,
    )
    baseline_bqs = conserved_densities_from_flavor(baseline_rho)
    masses = calculate_mass_vec(model, meanfield_state(baseline.x_state).phi)
    qp = (
        m=(u=Float64(masses[1]), d=Float64(masses[2]), s=Float64(masses[3])),
        μ=(u=Float64(baseline.mu_vec[1]), d=Float64(baseline.mu_vec[2]), s=Float64(baseline.mu_vec[3])),
        # Keep the equilibrium condensates beside the derived masses so the
        # legacy A/B post-processing can use the phi-native adapter.
        phi=(u=Float64(baseline.x_state[1]), d=Float64(baseline.x_state[2]), s=Float64(baseline.x_state[3])),
    )
    state = meanfield_state(baseline.x_state)
    tp = (T=T_fm, Φ=Float64(state.Phi), Φbar=Float64(state.PhiBar), ξ=0.0)
    meson_mu = charged_meson_chemical_potentials(conserved_mu.mu_Q, conserved_mu.mu_S)

    legacy = _charged_density_set(meson_mu, qp, tp, settings)
    kmt_params = _model_kmt_couplings(model)
    full_kernel = build_full_kmt_interaction(
        state.phi;
        G=kmt_params.G,
        K=kmt_params.K,
    )
    full = _charged_density_set(meson_mu, qp, tp, settings; interaction=full_kernel)

    statuses = (
        legacy.pi_plus.status, legacy.pi_minus.status, legacy.K_plus.status, legacy.K_minus.status,
        full.pi_plus.status, full.pi_minus.status, full.K_plus.status, full.K_minus.status,
    )
    converged = all(status === :ok for status in statuses)
    reason = converged ? :ok : :bu_density_status_not_ok
    message = converged ? "" : "one or more BU density channels returned a non-ok status"
    safe_ratio(numerator, denominator) = denominator > 0.0 ? numerator / denominator : NaN
    kernel_couplings = (
        K12_P=charged_coupling(full_kernel, :K12, :P),
        K45_P=charged_coupling(full_kernel, :K45, :P),
        K67_P=charged_coupling(full_kernel, :K67, :P),
        K03_P=full_kernel.neutral_P[1, 2],
        K38_P=full_kernel.neutral_P[2, 3],
        G=kmt_params.G,
        K=kmt_params.K,
        source=kmt_params.source,
    )
    legacy_payload = (
        n_pi_plus=Float64(legacy.pi_plus.density),
        n_pi_minus=Float64(legacy.pi_minus.density),
        n_K_plus=Float64(legacy.K_plus.density),
        n_K_minus=Float64(legacy.K_minus.density),
        Kplus_over_piplus=safe_ratio(legacy.K_plus.density, legacy.pi_plus.density),
        Kminus_over_piminus=safe_ratio(legacy.K_minus.density, legacy.pi_minus.density),
        pi_plus=legacy.pi_plus,
        pi_minus=legacy.pi_minus,
        K_plus=legacy.K_plus,
        K_minus=legacy.K_minus,
    )
    full_payload = (
        n_pi_plus=Float64(full.pi_plus.density),
        n_pi_minus=Float64(full.pi_minus.density),
        n_K_plus=Float64(full.K_plus.density),
        n_K_minus=Float64(full.K_minus.density),
        Kplus_over_piplus=safe_ratio(full.K_plus.density, full.pi_plus.density),
        Kminus_over_piminus=safe_ratio(full.K_minus.density, full.pi_minus.density),
        pi_plus=full.pi_plus,
        pi_minus=full.pi_minus,
        K_plus=full.K_plus,
        K_minus=full.K_minus,
    )

    return (
        converged=converged,
        reason=reason,
        message=message,
        baseline=baseline,
        baseline_mu=conserved_mu,
        baseline_bqs=baseline_bqs,
        meson_mu=meson_mu,
        qp=qp,
        thermo=tp,
        masses=(u=qp.m.u, d=qp.m.d, s=qp.m.s),
        legacy=legacy_payload,
        full=full_payload,
        full_kernel=full_kernel,
        kernel_couplings=kernel_couplings,
        target_ratio=target_ratio,
        rho_S_target=rho_S_target,
        rho0=rho0,
        baseline_elapsed_s=baseline_elapsed_s,
        legacy_density_elapsed_s=legacy.elapsed_s,
        full_density_elapsed_s=full.elapsed_s,
        total_elapsed_s=baseline_elapsed_s + legacy.elapsed_s + full.elapsed_s,
    )
end

end # module MesonConservedChargeFeedbackRuntime
