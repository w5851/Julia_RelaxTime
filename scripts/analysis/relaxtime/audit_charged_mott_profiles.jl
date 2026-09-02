"""
Pair real charged profiles on the two sides of a candidate Mott crossing.

This is an analysis-only wrapper around `audit_charged_phase_backend.jl`. It
does not infer a Mott temperature, promote a phase convention, or alter any
production route. The default temperatures are explicit diagnostic choices
and must be replaced by a resolved mass-threshold crossing before use as a
physical Mott result.
"""

const MOTT_PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "audit_charged_phase_backend.jl"))

using CSV
using .Models
using Main.RelaxTime.AFieldBuilder: build_A_triplet
using Main.RelaxTime.ChargedRPAKernel: charged_rpa_spec, charged_rpa_coupling
using Main.RelaxTime.ChargedRPAProvider: charged_pair_continuum_thresholds
using Main.RelaxTime.ChargedPhaseBackend: strict_mott_gate
using Main.RelaxTime.MesonInteractionKernel: build_full_kmt_interaction

const DEFAULT_OUTPUT = joinpath(
    MOTT_PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "analysis",
    "charged_rpa_phase_backend", "strict_mott_profile_pair.csv",
)

@inline _env_float_mott(name::AbstractString, default::Real) =
    parse(Float64, get(ENV, name, string(default)))
@inline _env_int_mott(name::AbstractString, default::Integer) =
    parse(Int, get(ENV, name, string(default)))

function _background_at(T_MeV::Real)
    ENV["CHARGED_PHASE_T_MEV"] = string(T_MeV)
    return _solve_background()
end

function _profile_at(background, meson::Symbol, settings)
    state = Models.meanfield_state(background.result.x_state)
    masses = (u=Float64(background.result.masses[1]), d=Float64(background.result.masses[2]),
              s=Float64(background.result.masses[3]))
    chemical_potentials = (u=Float64(background.result.mu_vec[1]),
                           d=Float64(background.result.mu_vec[2]),
                           s=Float64(background.result.mu_vec[3]))
    thermo = (T=background.T_fm, Φ=Float64(state.Phi), Φbar=Float64(state.PhiBar), ξ=0.0)
    A_values = build_A_triplet(
        (m=masses, μ=chemical_potentials), thermo;
        p_nodes=_env_int_mott("CHARGED_PHASE_A_NODES", 64),
        p_max=_env_float_mott("CHARGED_PHASE_A_PMAX", 16.0),
        use_aniso=false,
    )
    kernel = build_full_kmt_interaction(
        state.phi;
        G=background.model.params.G_fm2,
        K=background.model.params.K_fm5,
    )
    mass = meson === :pi_plus || meson === :pi_minus ?
        background.point.meson_results[:pi].mass : background.point.meson_results[:K].mass
    spec = charged_rpa_spec(meson)
    pair_threshold = charged_pair_continuum_thresholds(
        0.0,
        getproperty(masses, spec.pair[1]),
        getproperty(masses, spec.pair[2]),
        getproperty(chemical_potentials, spec.pair[1]),
        getproperty(chemical_potentials, spec.pair[2]),
    )
    threshold_q0 = Float64(pair_threshold.k0_threshold_inv_fm)
    lambda_threshold_q0 = Float64(pair_threshold.lambda_threshold_inv_fm)
    # A high-T pion threshold can fall below the ordinary phase audit window.
    # Lower only this diagnostic window so the threshold endpoint is sampled;
    # record the resulting setting rather than silently extrapolating it.
    omega_min = min(Float64(settings.omega_min), max(1.0e-3, 0.5 * threshold_q0))
    profile_settings = merge(settings, (omega_min=omega_min,))
    item = _run_mode(
        meson, mass, masses, chemical_potentials, thermo, A_values, kernel, profile_settings,
    )
    isempty(item.result.q_profiles) && error("no q profiles returned for $(meson)")
    first_profile = item.result.q_profiles[1]
    first_profile.gate === nothing && error("Mott audit requires a phase gate")
    return (
        item=item,
        profile=first_profile.profile,
        gate=first_profile.gate,
        q=Float64(first_profile.q),
        threshold=Float64(first_profile.threshold),
        bound_state_count=Int(first_profile.bound_state_count),
        mass=Float64(mass),
        threshold_q0=threshold_q0,
        lambda_threshold_q0=lambda_threshold_q0,
        bound_state_q_count_min=minimum(Int(profile.bound_state_count) for profile in item.result.q_profiles),
        bound_state_q_count_max=maximum(Int(profile.bound_state_count) for profile in item.result.q_profiles),
        bound_state_complex_q_count=count(
            profile -> profile.bound_state_diagnostic !== nothing &&
                       profile.bound_state_diagnostic.status === :complex_subthreshold,
            item.result.q_profiles,
        ),
        bound_state_q_points=length(item.result.q_profiles),
    )
end

function main()
    T_before = _env_float_mott("CHARGED_PHASE_MOTT_T_BEFORE_MEV", 170.0)
    T_after = _env_float_mott("CHARGED_PHASE_MOTT_T_AFTER_MEV", 230.0)
    settings = _settings(false)
    rows = NamedTuple[]
    for meson in CHARGED_MODES
        before = _profile_at(_background_at(T_before), meson, settings)
        after = _profile_at(_background_at(T_after), meson, settings)
        transition = strict_mott_gate(
            before.profile, after.profile;
            before_threshold=before.threshold,
            after_threshold=after.threshold,
            before_bound_state_count=before.bound_state_count,
            after_bound_state_count=after.bound_state_count,
            expected_bound_state_drop=_env_int_mott("CHARGED_PHASE_MOTT_EXPECTED_DROP", 1),
        )
        before_diag = before.item.result.q_profiles[1].bound_state_diagnostic
        after_diag = after.item.result.q_profiles[1].bound_state_diagnostic
        push!(rows, (
            meson=String(meson), T_before_MeV=T_before, T_after_MeV=T_after,
            q_inv_fm=before.q,
            before_mass_inv_fm=before.mass, after_mass_inv_fm=after.mass,
            before_threshold_q0_inv_fm=before.threshold_q0,
            after_threshold_q0_inv_fm=after.threshold_q0,
            before_pair_threshold_lambda_q0_inv_fm=before.lambda_threshold_q0,
            after_pair_threshold_lambda_q0_inv_fm=after.lambda_threshold_q0,
            before_mott_gap_inv_fm=before.mass - before.threshold_q0,
            after_mott_gap_inv_fm=after.mass - after.threshold_q0,
            before_bound_state_count=before.bound_state_count,
            after_bound_state_count=after.bound_state_count,
            before_bound_state_status=String(before_diag.status),
            after_bound_state_status=String(after_diag.status),
            before_bound_state_q_count_min=before.bound_state_q_count_min,
            before_bound_state_q_count_max=before.bound_state_q_count_max,
            after_bound_state_q_count_min=after.bound_state_q_count_min,
            after_bound_state_q_count_max=after.bound_state_q_count_max,
            before_bound_state_complex_q_count=before.bound_state_complex_q_count,
            after_bound_state_complex_q_count=after.bound_state_complex_q_count,
            bound_state_q_points=before.bound_state_q_points,
            before_levinson_passed=Bool(before.gate.levinson.passed),
            after_levinson_passed=Bool(after.gate.levinson.passed),
            before_phase_gate_passed=Bool(before.gate.passed),
            after_phase_gate_passed=Bool(after.gate.passed),
            mott_gate_passed=Bool(transition.passed),
            transition_phase_residual=Float64(transition.transition.phase_residual),
            expected_bound_state_drop=Int(transition.transition.expected_bound_state_drop),
            production_candidate_status="not_authorized",
        ))
    end
    output = get(ENV, "CHARGED_PHASE_MOTT_OUTPUT", DEFAULT_OUTPUT)
    mkpath(dirname(output))
    CSV.write(output, rows)
    println("[charged-mott-profiles] wrote $(length(rows)) rows to $(output)")
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
