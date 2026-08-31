"""
Diagnostic fixed-background comparison for the charged-RPA provider.

The script solves one `FixedMuBConservedCharges` quark-only background and
compares the strict ordered retarded bubble with both legacy real-axis
oracles. It does not change any production entrypoint and does not perform a
Beth--Uhlenbeck density integral.
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using CSV
using Printf: @printf
using .Models
using Main.Constants_PNJL: ħc_MeV_fm
using Main.RelaxTime.AFieldBuilder: build_A_triplet
using Main.RelaxTime.ChargedRPAKernel: charged_rpa_spec, charged_rpa_coupling,
                                             charged_rpa_inverse, charged_rpa_propagator
using Main.RelaxTime.ChargedRPAProvider: charged_polarization
using Main.RelaxTime.MesonInteractionKernel: build_full_kmt_interaction

const FIXEDPOINT_MESONS = (:pi_plus, :pi_minus, :K_plus, :K_minus)
const DEFAULT_OUTPUT = joinpath(
    PROJECT_ROOT,
    "data", "outputs", "results", "relaxtime", "analysis",
    "charged_rpa_ordered_fixedpoint", "t170_mub0_mub240.csv",
)

@inline _env_float(name::AbstractString, default::Real) = parse(Float64, get(ENV, name, string(default)))
@inline _env_int(name::AbstractString, default::Integer) = parse(Int, get(ENV, name, string(default)))

function _route_settings()
    coarse_eta = _env_float("CHARGED_RPA_COARSE_ETA", 0.01)
    refined_eta = _env_float("CHARGED_RPA_REFINED_ETA", 0.005)
    coarse_nodes = _env_int("CHARGED_RPA_COARSE_ENERGY_NODES", 512)
    refined_nodes = _env_int("CHARGED_RPA_REFINED_ENERGY_NODES", 1024)
    return (
        (label=:ordered_retarded_eta_coarse, prescription=:ordered_retarded, eta=coarse_eta, nodes=refined_nodes),
        (label=:ordered_retarded_nodes_coarse, prescription=:ordered_retarded, eta=refined_eta, nodes=coarse_nodes),
        (label=:ordered_retarded_refined, prescription=:ordered_retarded, eta=refined_eta, nodes=refined_nodes),
        (label=:ordered_legacy_B0, prescription=:ordered_legacy_B0, eta=0.0, nodes=4),
        (label=:legacy_symmetrized_B0, prescription=:legacy_symmetrized_B0, eta=0.0, nodes=4),
    )
end

function _background_specs()
    T_MeV = _env_float("CHARGED_RPA_FIXED_T_MEV", 170.0)
    return (
        (label=:charge_symmetric_reference, T_MeV=T_MeV, muB_MeV=0.0),
        (
            label=:finite_bqs_reference,
            T_MeV=T_MeV,
            muB_MeV=_env_float("CHARGED_RPA_FIXED_MUB_MEV", 240.0),
        ),
    )
end

@inline function _relative_complex_difference(value::Complex, reference::Complex)
    scale = max(abs(reference), eps(Float64))
    return abs(value - reference) / scale
end

function _probe_specs(meson::Symbol, masses)
    spec = charged_rpa_spec(meson)
    m1 = getproperty(masses, spec.pair[1])
    m2 = getproperty(masses, spec.pair[2])
    threshold = m1 + m2
    return (
        (kind=:below_threshold_q0, omega=0.85 * threshold, q=0.0, threshold=threshold),
        (kind=:above_threshold_finite_q, omega=sqrt(0.35^2 + (threshold + 0.30)^2), q=0.35, threshold=threshold),
    )
end

function _solve_background(spec)
    T_MeV = Float64(spec.T_MeV)
    muB_MeV = Float64(spec.muB_MeV)
    p_num = _env_int("CHARGED_RPA_FIXED_P_NUM", 8)
    t_num = _env_int("CHARGED_RPA_FIXED_T_NUM", 4)
    residual_norm_max = _env_float("CHARGED_RPA_FIXED_RESIDUAL_MAX", 1e-5)
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
        residual_norm_max=residual_norm_max,
        iterations=200,
        xi=0.0,
    )
    result.converged || error(
        "fixed BQS background did not converge: residual=$(result.residual_norm)",
    )
    return (
        model=model,
        result=result,
        label=spec.label,
        T_MeV=T_MeV,
        muB_MeV=muB_MeV,
        T_fm=T_fm,
        muB_fm=muB_fm,
        p_num=p_num,
        t_num=t_num,
        residual_norm_max=residual_norm_max,
    )
end

@inline function _charge_partner(meson::Symbol)
    meson === :pi_plus && return :pi_minus
    meson === :pi_minus && return :pi_plus
    meson === :K_plus && return :K_minus
    meson === :K_minus && return :K_plus
    throw(ArgumentError("unsupported charged meson $(meson)"))
end

function _fixedpoint_rows(background)
    result = background.result
    state = Models.meanfield_state(result.x_state)
    masses = (
        u=Float64(result.masses[1]),
        d=Float64(result.masses[2]),
        s=Float64(result.masses[3]),
    )
    chemical_potentials = (
        u=Float64(result.mu_vec[1]),
        d=Float64(result.mu_vec[2]),
        s=Float64(result.mu_vec[3]),
    )
    thermo = (
        T=Float64(background.T_fm),
        Φ=Float64(state.Phi),
        Φbar=Float64(state.PhiBar),
        ξ=0.0,
    )
    quark_params = (m=masses, μ=chemical_potentials)
    A_nodes = _env_int("CHARGED_RPA_A_NODES", 128)
    A_pmax = _env_float("CHARGED_RPA_A_PMAX", 20.0)
    A_values = build_A_triplet(
        quark_params,
        thermo;
        p_nodes=A_nodes,
        p_max=A_pmax,
        use_aniso=false,
    )
    kernel = build_full_kmt_interaction(
        state.phi;
        G=background.model.params.G_fm2,
        K=background.model.params.K_fm5,
    )
    conserved_mu = Models.conserved_mu_from_flavor(result.mu_vec...)
    rows = NamedTuple[]

    for meson in FIXEDPOINT_MESONS
        spec = charged_rpa_spec(meson)
        coupling = charged_rpa_coupling(kernel, spec)
        for probe in _probe_specs(meson, masses)
            evaluations = NamedTuple[]
            for route in _route_settings()
                evaluated = charged_polarization(
                    spec,
                    probe.omega,
                    probe.q,
                    masses,
                    chemical_potentials,
                    thermo,
                    A_values;
                    prescription=route.prescription,
                    eta_inv_fm=route.eta,
                    energy_nodes=route.nodes,
                )
                inverse = charged_rpa_inverse(spec, coupling, evaluated.value)
                propagator = charged_rpa_propagator(spec, coupling, evaluated.value)
                push!(evaluations, (
                    route=route,
                    evaluated=evaluated,
                    inverse=ComplexF64(inverse),
                    propagator=ComplexF64(propagator),
                ))
            end
            reference = only(filter(
                item -> item.route.label === :ordered_retarded_refined,
                evaluations,
            ))

            for item in evaluations
                evaluated = item.evaluated
                push!(rows, (
                    route="charged_rpa_ordered_fixedpoint_diagnostic",
                    production_candidate_status="not_authorized",
                    background_label=String(background.label),
                    T_MeV=background.T_MeV,
                    muB_MeV=background.muB_MeV,
                    target_Q_over_B=0.4,
                    charge_ratio_constraint_status=abs(background.muB_MeV) <= 1e-12 ?
                        "zero_density_degenerate_reference" : "active",
                    rho_S_target_fm3=0.0,
                    gap_converged=result.converged,
                    gap_residual_norm=Float64(result.residual_norm),
                    p_num=background.p_num,
                    t_num=background.t_num,
                    muB_solved_MeV=Float64(conserved_mu.mu_B * ħc_MeV_fm),
                    muQ_solved_MeV=Float64(conserved_mu.mu_Q * ħc_MeV_fm),
                    muS_solved_MeV=Float64(conserved_mu.mu_S * ħc_MeV_fm),
                    phi_u_fm3=Float64(state.phi[1]),
                    phi_d_fm3=Float64(state.phi[2]),
                    phi_s_fm3=Float64(state.phi[3]),
                    Phi=thermo.Φ,
                    PhiBar=thermo.Φbar,
                    m_u_inv_fm=masses.u,
                    m_d_inv_fm=masses.d,
                    m_s_inv_fm=masses.s,
                    mu_u_inv_fm=chemical_potentials.u,
                    mu_d_inv_fm=chemical_potentials.d,
                    mu_s_inv_fm=chemical_potentials.s,
                    A_u_inv_fm2=Float64(A_values.u),
                    A_d_inv_fm2=Float64(A_values.d),
                    A_s_inv_fm2=Float64(A_values.s),
                    A_nodes=A_nodes,
                    A_pmax_inv_fm=A_pmax,
                    meson=String(meson),
                    flavor1=String(spec.pair[1]),
                    flavor2=String(spec.pair[2]),
                    kernel_pair=String(spec.kernel_pair),
                    coupling_fm2=Float64(coupling),
                    probe_kind=String(probe.kind),
                    omega_inv_fm=Float64(probe.omega),
                    q_inv_fm=Float64(probe.q),
                    threshold_inv_fm=Float64(probe.threshold),
                    variant=String(item.route.label),
                    prescription=String(evaluated.prescription),
                    provider=String(evaluated.provider),
                    analytic_scope=String(evaluated.analytic_scope),
                    num_s_quark=evaluated.num_s_quark,
                    eta_inv_fm=Float64(evaluated.eta_inv_fm),
                    energy_nodes=Int(evaluated.energy_nodes),
                    Pi_real_inv_fm2=Float64(real(evaluated.value)),
                    Pi_imag_inv_fm2=Float64(imag(evaluated.value)),
                    inverse_real=Float64(real(item.inverse)),
                    inverse_imag=Float64(imag(item.inverse)),
                    propagator_real_fm2=Float64(real(item.propagator)),
                    propagator_imag_fm2=Float64(imag(item.propagator)),
                    Pi_relative_to_refined=_relative_complex_difference(
                        evaluated.value,
                        reference.evaluated.value,
                    ),
                    inverse_relative_to_refined=_relative_complex_difference(
                        item.inverse,
                        reference.inverse,
                    ),
                    propagator_relative_to_refined=_relative_complex_difference(
                        item.propagator,
                        reference.propagator,
                    ),
                ))
            end
        end
    end

    return map(rows) do row
        partner_name = String(_charge_partner(Symbol(row.meson)))
        partner = only(filter(rows) do candidate
            candidate.meson == partner_name &&
                candidate.probe_kind == row.probe_kind &&
                candidate.variant == row.variant
        end)
        row_Pi = ComplexF64(row.Pi_real_inv_fm2, row.Pi_imag_inv_fm2)
        partner_Pi = ComplexF64(partner.Pi_real_inv_fm2, partner.Pi_imag_inv_fm2)
        row_propagator = ComplexF64(row.propagator_real_fm2, row.propagator_imag_fm2)
        partner_propagator = ComplexF64(
            partner.propagator_real_fm2,
            partner.propagator_imag_fm2,
        )
        charge_partner_gate = abs(background.muB_MeV) <= 1e-12
        return merge(row, (
            positive_energy_charge_partner=partner_name,
            positive_energy_charge_partner_status=charge_partner_gate ?
                "evaluated_muB0" : "not_applicable_finite_mu",
            Pi_positive_energy_partner_relative=charge_partner_gate ?
                _relative_complex_difference(row_Pi, partner_Pi) : NaN,
            propagator_positive_energy_partner_relative=charge_partner_gate ?
                _relative_complex_difference(row_propagator, partner_propagator) : NaN,
        ))
    end
end

function build_fixedpoint_rows()
    rows = NamedTuple[]
    for spec in _background_specs()
        append!(rows, _fixedpoint_rows(_solve_background(spec)))
    end
    return rows
end

function main()
    output_path = get(ENV, "CHARGED_RPA_FIXEDPOINT_OUTPUT", DEFAULT_OUTPUT)
    rows = build_fixedpoint_rows()
    mkpath(dirname(output_path))
    CSV.write(output_path, rows)
    @printf("wrote %d charged-RPA fixed-point rows to %s\n", length(rows), output_path)
    return output_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
