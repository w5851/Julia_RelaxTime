"""
Physical fixed-background diagnostic for the strict charged phase/BU backend.

The script evaluates the four charge-resolved modes on one finite-BQS
quark-only background using the ordered retarded `ChargedRPAProvider` bubble
and the full KMT charged coupling.  It is deliberately analysis-only: failed
endpoint, Levinson, Bose-support, or density-sign gates are written as rows
and never converted to production values.
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using CSV
using .Models
using Main.Constants_PNJL: ħc_MeV_fm
using Main.RelaxTime.AFieldBuilder: build_A_triplet
using Main.RelaxTime.ChargedRPAKernel: charged_rpa_spec, charged_rpa_coupling
using Main.RelaxTime.ChargedRPAProvider: charged_polarization
using Main.RelaxTime.ChargedRPAKernel: charged_rpa_inverse
using Main.RelaxTime.BUPhaseGates: count_bound_states
using Main.RelaxTime.ChargedPhaseBackend: StrictChargedPhaseSpec,
                                           strict_charged_rpa_bu_density,
                                           strict_density_convergence_gate
using Main.RelaxTime.MesonInteractionKernel: build_full_kmt_interaction

const CHARGED_MODES = (:pi_plus, :pi_minus, :K_plus, :K_minus)
const DEFAULT_OUTPUT = joinpath(
    PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "analysis",
    "charged_rpa_phase_backend", "strict_fixed_bqs_t170_mub240.csv",
)

@inline _env_float(name::AbstractString, default::Real) =
    parse(Float64, get(ENV, name, string(default)))
@inline _env_int(name::AbstractString, default::Integer) =
    parse(Int, get(ENV, name, string(default)))
# Supported diagnostic choices: :ordered_retarded and :ordered_pv_cut.
@inline _prescription() = Symbol(get(ENV, "CHARGED_PHASE_PRESCRIPTION", "ordered_retarded"))
@inline _omega_measure() = Symbol(get(ENV, "CHARGED_PHASE_OMEGA_MEASURE", "single_charge_domega_over_pi"))

function _solve_background()
    T_MeV = _env_float("CHARGED_PHASE_T_MEV", 170.0)
    muB_MeV = _env_float("CHARGED_PHASE_MUB_MEV", 240.0)
    p_num = _env_int("CHARGED_PHASE_P_NUM", 8)
    t_num = _env_int("CHARGED_PHASE_T_NUM", 4)
    T_fm = T_MeV / ħc_MeV_fm
    muB_fm = muB_MeV / ħc_MeV_fm
    model = Models.create_model(:PNJL)
    mode = Models.FixedMuBConservedCharges(muB_fm, 0.4, 0.0)
    result = Models.solve(
        model,
        mode,
        T_fm;
        p_num=p_num,
        t_num=t_num,
        residual_norm_max=_env_float("CHARGED_PHASE_GAP_RESIDUAL_MAX", 1.0e-5),
        iterations=_env_int("CHARGED_PHASE_GAP_ITERATIONS", 200),
        xi=0.0,
    )
    result.converged || error("finite-BQS background did not converge: residual=$(result.residual_norm)")
    point = Models.solve_meson_point_from_equilibrium(
        result,
        T_fm;
        mesons=(:pi, :K),
        mass_kwargs=(;iterations=20),
        mixed_branch_align=:strict_sign_binding,
    )
    return (model=model, result=result, point=point, T_MeV=T_MeV, muB_MeV=muB_MeV,
            T_fm=T_fm, p_num=p_num, t_num=t_num)
end

function _settings(refined::Bool)
    suffix = refined ? "_REFINED" : ""
    get_float(name, default) = _env_float("CHARGED_PHASE_" * name * suffix, default)
    get_int(name, default) = _env_int("CHARGED_PHASE_" * name * suffix, default)
    return (
        eta=get_float("ETA_INV_FM", refined ? 5.0e-3 : 1.0e-2),
        polarization_nodes=get_int("POLARIZATION_NODES", refined ? 32 : 16),
        qmax=get_float("QMAX", refined ? 2.0 : 1.5),
        q_nodes=get_int("Q_NODES", refined ? 4 : 2),
        omega_min=get_float("OMEGA_MIN", 0.5),
        omega_max=get_float("OMEGA_MAX", refined ? 10.0 : 8.0),
        omega_nodes=get_int("OMEGA_NODES", refined ? 24 : 12),
        bound_state_nodes=get_int("BOUND_STATE_NODES", refined ? 96 : 64),
        variant=refined ? :refined : :coarse,
    )
end

function _run_mode(meson::Symbol, meson_mass::Real, masses, chemical_potentials, thermo, A_values, kernel, settings)
    spec = charged_rpa_spec(meson)
    coupling = charged_rpa_coupling(kernel, spec)
    prescription = _prescription()
    m1 = getproperty(masses, spec.pair[1])
    m2 = getproperty(masses, spec.pair[2])
    threshold_fn = q -> sqrt(q^2 + m1^2) + sqrt(q^2 + m2^2)
    polarization_fn = (ω, q) -> charged_polarization(
        spec,
        ω,
        q,
        masses,
        chemical_potentials,
        thermo,
        A_values;
        prescription=prescription,
        eta_inv_fm=settings.eta,
        energy_nodes=settings.polarization_nodes,
    ).value
    inverse_fn = (ω, q) -> charged_rpa_inverse(
        spec,
        coupling,
        polarization_fn(ω, q),
    )
    phase_spec = StrictChargedPhaseSpec(
        tail_points=_env_int("CHARGED_PHASE_TAIL_POINTS", 4),
        tail_tolerance=_env_float("CHARGED_PHASE_TAIL_TOLERANCE", 0.2π),
    )
    μ_meson = meson === :pi_plus ? chemical_potentials.u - chemical_potentials.d :
        meson === :pi_minus ? chemical_potentials.d - chemical_potentials.u :
        meson === :K_plus ? chemical_potentials.u - chemical_potentials.s :
        chemical_potentials.s - chemical_potentials.u
    result = strict_charged_rpa_bu_density(
        spec,
        coupling,
        polarization_fn,
        meson_mass,
        thermo.T;
        μ=μ_meson,
        qmax=settings.qmax,
        q_nodes=settings.q_nodes,
        omega_min=settings.omega_min,
        omega_max=settings.omega_max,
        omega_nodes=settings.omega_nodes,
        threshold=threshold_fn,
        # Count subthreshold real-axis zeros independently.  A finite eta
        # makes this return :complex_subthreshold; that failure remains
        # visible in q_profiles and is not converted into a physical count.
        bound_state_count=q -> count_bound_states(
            inverse_fn,
            q,
            threshold_fn(q);
            omega_nodes=settings.bound_state_nodes,
            imag_tolerance=_env_float("CHARGED_PHASE_BOUND_STATE_IMAG_TOL", 1.0e-8),
        ),
        phase_spec=phase_spec,
        omega_measure=_omega_measure(),
        require_levinson=true,
    )
    return (meson=meson, spec=spec, coupling=coupling, μ_meson=μ_meson, result=result, settings=settings)
end

function main()
    background = _solve_background()
    state = Models.meanfield_state(background.result.x_state)
    masses = (u=Float64(background.result.masses[1]), d=Float64(background.result.masses[2]),
              s=Float64(background.result.masses[3]))
    chemical_potentials = (u=Float64(background.result.mu_vec[1]),
                           d=Float64(background.result.mu_vec[2]),
                           s=Float64(background.result.mu_vec[3]))
    thermo = (T=background.T_fm, Φ=Float64(state.Phi), Φbar=Float64(state.PhiBar), ξ=0.0)
    A_values = build_A_triplet(
        (m=masses, μ=chemical_potentials),
        thermo;
        p_nodes=_env_int("CHARGED_PHASE_A_NODES", 64),
        p_max=_env_float("CHARGED_PHASE_A_PMAX", 16.0),
        use_aniso=false,
    )
    kernel = build_full_kmt_interaction(
        state.phi;
        G=background.model.params.G_fm2,
        K=background.model.params.K_fm5,
    )
    rows = NamedTuple[]
    for meson in CHARGED_MODES
        meson_mass = meson === :pi_plus || meson === :pi_minus ?
            background.point.meson_results[:pi].mass : background.point.meson_results[:K].mass
        item = _run_mode(meson, meson_mass, masses, chemical_potentials, thermo, A_values, kernel, _settings(false))
        refined_item = _run_mode(meson, meson_mass, masses, chemical_potentials, thermo, A_values, kernel, _settings(true))
        convergence = strict_density_convergence_gate(item.result, refined_item.result; rtol=0.05)
        for evaluated_item in (item, refined_item)
            result = evaluated_item.result
            settings = evaluated_item.settings
            q_profiles = hasproperty(result, :q_profiles) ? result.q_profiles : NamedTuple[]
            tail_failed_q_count = count(profile -> !Bool(profile.profile.tail_stable), q_profiles)
            root_failed_q_count = count(
                profile -> profile.gate !== nothing && !Bool(profile.gate.roots.passed),
                q_profiles,
            )
            levinson_failed_q_count = count(
                profile -> profile.gate !== nothing && !Bool(profile.gate.levinson.passed),
                q_profiles,
            )
            bound_state_diagnostics = [
                profile.bound_state_diagnostic for profile in q_profiles
                if profile.bound_state_diagnostic !== nothing
            ]
            bound_state_complex_q_count = count(
                diagnostic -> diagnostic.status === :complex_subthreshold,
                bound_state_diagnostics,
            )
            bound_state_status = isempty(bound_state_diagnostics) ? "not_provided" :
                (bound_state_complex_q_count > 0 ? "complex_subthreshold" :
                 (all(diagnostic -> diagnostic.status === :ok, bound_state_diagnostics) ? "ok" : "diagnostic_failed"))
            gate_profiles = [profile.gate for profile in q_profiles if profile.gate !== nothing]
            levinson_residual_values = [
                abs(Float64(gate.levinson.levinson_residual)) for gate in gate_profiles
            ]
            levinson_residual_max = isempty(levinson_residual_values) ? NaN : maximum(levinson_residual_values)
            threshold_phase_values = [Float64(gate.levinson.threshold_phase) for gate in gate_profiles]
            threshold_phase_min = isempty(threshold_phase_values) ? NaN : minimum(threshold_phase_values)
            threshold_phase_max = isempty(threshold_phase_values) ? NaN : maximum(threshold_phase_values)
            tail_span_values = [Float64(gate.levinson.tail_span) for gate in gate_profiles]
            tail_span_max = isempty(tail_span_values) ? NaN : maximum(tail_span_values)
            push!(rows, (
            route="strict_charged_phase_backend_fixed_bqs",
            production_candidate_status="not_authorized",
            meson=String(meson),
            variant=String(settings.variant),
            mu_meson_inv_fm=Float64(evaluated_item.μ_meson),
            pair="$(evaluated_item.spec.pair[1])bar$(evaluated_item.spec.pair[2])",
            kernel_pair=String(evaluated_item.spec.kernel_pair),
            polarization_prescription=String(_prescription()),
            T_MeV=background.T_MeV,
            muB_MeV=background.muB_MeV,
            gap_converged=Bool(background.result.converged),
            gap_residual_norm=Float64(background.result.residual_norm),
            p_num=background.p_num,
            t_num=background.t_num,
            meson_mass_inv_fm=Float64(meson === :pi_plus || meson === :pi_minus ?
                background.point.meson_results[:pi].mass : background.point.meson_results[:K].mass),
            continuum_threshold_q0_inv_fm=Float64(getproperty(masses, evaluated_item.spec.pair[1]) + getproperty(masses, evaluated_item.spec.pair[2])),
            coupling_fm2=Float64(evaluated_item.coupling),
            density_fm3=Float64(result.density),
            accepted=Bool(result.accepted),
            status=String(result.status),
            omega_measure=String(result.omega_measure),
            omega_measure_factor=Float64(result.omega_measure_factor),
            eta_inv_fm=Float64(settings.eta),
            polarization_nodes=Int(settings.polarization_nodes),
            qmax_inv_fm=Float64(result.qmax),
            q_nodes=Int(result.q_nodes),
            omega_min_inv_fm=Float64(result.omega_min),
            omega_max_inv_fm=Float64(result.omega_max),
            omega_nodes=Int(result.omega_nodes),
            bound_state_nodes=Int(settings.bound_state_nodes),
            failed_q_count=Int(result.failed_q_count),
            tail_failed_q_count=tail_failed_q_count,
            root_failed_q_count=root_failed_q_count,
            levinson_failed_q_count=levinson_failed_q_count,
            bound_state_complex_q_count=bound_state_complex_q_count,
            bound_state_status=bound_state_status,
            levinson_residual_max=levinson_residual_max,
            threshold_phase_min=threshold_phase_min,
            threshold_phase_max=threshold_phase_max,
            tail_span_max=tail_span_max,
            density_finite=Bool(result.density_finite),
            density_nonnegative=Bool(result.density_nonnegative),
            threshold_policy=String(result.threshold_policy),
            bound_state_policy=String(result.bound_state_policy),
            convergence_passed=Bool(convergence.passed),
            convergence_relative_difference=convergence.numeric === nothing ? NaN : Float64(convergence.numeric.relative_difference),
            coarse_density_fm3=Float64(item.result.density),
            refined_density_fm3=Float64(refined_item.result.density),
            ))
        end
    end
    output = get(ENV, "CHARGED_PHASE_OUTPUT", DEFAULT_OUTPUT)
    mkpath(dirname(output))
    CSV.write(output, rows)
    println("[charged-phase-backend] wrote $(length(rows)) rows to $(output)")
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
