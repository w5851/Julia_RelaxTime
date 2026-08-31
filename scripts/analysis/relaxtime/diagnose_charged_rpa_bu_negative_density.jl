"""
Fixed-background diagnostic for the charged-RPA/BU negative-density audit.

This script intentionally stays outside the production workflow.  It solves one
finite-BQS quark-only point, obtains the existing pi/K mass point, and records:

* four density algorithms on the same state;
* q-pole acceptance arrays;
* raw/anchored `arg(D)` and `-arg(D^{-1})` phase profiles;
* the effect of `current`/`gbu_reference`, phase convention, anchor, and
  positive-energy omega measure.

The output is a local diagnostic text file.  No production default or baseline
is changed.
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "workflow_apps", "WorkflowParamAdapters.jl"))

using .Models
using .WorkflowParamAdapters: normalize_quark_params, normalize_thermo_params
using Printf: @sprintf
using Main.Constants_PNJL: ħc_MeV_fm
using Main.RelaxTime.GaussLegendre: gauleg
using Main.RelaxTime.AFieldBuilder: ensure_quark_params_has_A
using Main.RelaxTime.BUPhaseGates: anchor_phase_high_energy
using Main.RelaxTime.BUPhaseGates: bose_support_gate
using Main.RelaxTime.MesonDensity: phase_shift_meson_number_density,
                                    strict_bw_qpole_meson_number_density,
                                    stable_meson_number_density,
                                    strict_bw_meson_number_density,
                                    bose_distribution,
                                    _build_k_coeffs,
                                    _propagator_components,
                                    _propagator_phase,
                                    _phase_shift_weighted_phase

const OUTPUT_PATH = get(ENV, "CHARGED_RPA_NEGATIVE_DENSITY_OUTPUT", joinpath(
    PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "analysis",
    "charged_rpa_bu_negative_density", "t170_mub240_diagnostic.txt",
))

@inline _f(x) = @sprintf("%.12g", Float64(x))
@inline _complex_f(z) = string("(", _f(real(z)), ",", _f(imag(z)), ")")

function _solve_state()
    T_fm = 170.0 / ħc_MeV_fm
    muB_fm = 240.0 / ħc_MeV_fm
    model = Models.create_model(:PNJL)
    eq = Models.solve(
        model,
        Models.FixedMuBConservedCharges(muB_fm, 0.4, 0.0),
        T_fm;
        p_num=8,
        t_num=4,
        residual_norm_max=1e-5,
        iterations=200,
        xi=0.0,
    )
    eq.converged || error("finite-BQS background failed: residual=$(eq.residual_norm)")
    point = Models.solve_meson_point_from_equilibrium(
        eq,
        T_fm;
        mesons=(:pi, :K),
        mass_kwargs=(;iterations=20),
        mixed_branch_align=:strict_sign_binding,
    )
    qp = ensure_quark_params_has_A(
        normalize_quark_params(point.quark_params),
        normalize_thermo_params(point.thermo_params),
        p_nodes=128,
        p_max=20.0,
        use_aniso=false,
    )
    tp = normalize_thermo_params(point.thermo_params)
    return model, eq, point, qp, tp
end

function _mu_pairs(qp)
    return (
        pi_plus=Float64(qp.μ.u - qp.μ.d),
        pi_minus=Float64(qp.μ.d - qp.μ.u),
        K_plus=Float64(qp.μ.u - qp.μ.s),
        K_minus=Float64(qp.μ.s - qp.μ.u),
    )
end

function _density_suite(point, qp, tp, io)
    println(io, "[density_algorithms]")
    # Use omega_min=0.4 for this finite-BQS comparison so every charged Bose
    # window is normal-domain safe, including K+ with positive chemical potential.
    μs = _mu_pairs(qp)
    for charge in (:plus, :minus)
        pi_channel = charge === :plus ? :pi_plus : :pi_minus
        k_channel = charge === :plus ? :K_plus : :K_minus
        μ_pi = getproperty(μs, pi_channel)
        μ_K = getproperty(μs, k_channel)
        stable = stable_meson_number_density(
            Float64(point.meson_results[:pi].mass), tp.T;
            μ=μ_pi, degeneracy=1, qmax=4.0, num_q_nodes=64,
        )
        stable_k = stable_meson_number_density(
            Float64(point.meson_results[:K].mass), tp.T;
            μ=μ_K, degeneracy=1, qmax=4.0, num_q_nodes=64,
        )
        reduced = strict_bw_meson_number_density(
            Float64(point.meson_results[:pi].mass), Float64(point.meson_results[:pi].gamma), tp.T;
            μ=μ_pi, degeneracy=1, qmax=4.0, q_nodes=4, omega_min=0.4,
            omega_max=3.0, omega_nodes=8,
        )
        reduced_k = strict_bw_meson_number_density(
            Float64(point.meson_results[:K].mass), Float64(point.meson_results[:K].gamma), tp.T;
            μ=μ_K, degeneracy=1, qmax=4.0, q_nodes=4, omega_min=0.4,
            omega_max=3.0, omega_nodes=8,
        )
        qpole = try
            strict_bw_qpole_meson_number_density(
                pi_channel,
                Float64(point.meson_results[:pi].mass), Float64(point.meson_results[:pi].gamma),
                qp, tp; μ=μ_pi, degeneracy=1, qmax=4.0, q_nodes=4, omega_max=3.0,
                omega_nodes=8, solver_iterations=8, pole_residual_norm_max=1e-5,
                pole_require_converged=false,
            )
        catch err
            (density=NaN, accepted_flags=Bool[], converged_flags=Bool[], residual_norms=Float64[],
             error=sprint(showerror, err))
        end
        qpole_k = try
            strict_bw_qpole_meson_number_density(
                k_channel,
                Float64(point.meson_results[:K].mass), Float64(point.meson_results[:K].gamma),
                qp, tp; μ=μ_K, degeneracy=1, qmax=4.0, q_nodes=4, omega_max=3.0,
                omega_nodes=8, solver_iterations=8, pole_residual_norm_max=1e-5,
                pole_require_converged=false,
            )
        catch err
            (density=NaN, accepted_flags=Bool[], converged_flags=Bool[], residual_norms=Float64[],
             error=sprint(showerror, err))
        end
        println(io, "charge=$(charge) mu_pi=$(_f(μ_pi)) mu_K=$(_f(μ_K))")
        pi_support = bose_support_gate(Float64(point.meson_results[:pi].mass), μ_pi; omega_min=0.0, omega_max=3.0)
        k_support = bose_support_gate(Float64(point.meson_results[:K].mass), μ_K; omega_min=0.0, omega_max=3.0)
        println(io, "qpole_bose_support pi_status=$(pi_support.status) K_status=$(k_support.status) " *
            "pi_min_E_minus_mu=$(_f(pi_support.min_E_minus_mu)) K_min_E_minus_mu=$(_f(k_support.min_E_minus_mu)) " *
            "pi_window_safe=$(pi_support.integration_window_safe) K_window_safe=$(k_support.integration_window_safe)")
        println(io, "stable n_pi=$(_f(stable)) n_K=$(_f(stable_k)) ratio=$(_f(stable_k/stable))")
        println(io, "reduced_bw n_pi=$(_f(reduced.density)) n_K=$(_f(reduced_k.density)) ratio=$(_f(reduced_k.density/reduced.density))")
        println(io, "qpole n_pi=$(_f(qpole.density)) n_K=$(_f(qpole_k.density)) ratio=$(_f(qpole_k.density/qpole.density))")
        println(io, "qpole_pi accepted=$(qpole.accepted_flags) converged=$(qpole.converged_flags) residuals=$(qpole.residual_norms) error=$(get(qpole, :error, ""))")
        println(io, "qpole_K accepted=$(qpole_k.accepted_flags) converged=$(qpole_k.converged_flags) residuals=$(qpole_k.residual_norms) error=$(get(qpole_k, :error, ""))")
        println(io, "qpole_pi q=$(get(qpole, :q_values, Float64[])) E=$(get(qpole, :E_values, Float64[])) gamma=$(get(qpole, :gamma_values, Float64[]))")
        println(io, "qpole_K q=$(get(qpole_k, :q_values, Float64[])) E=$(get(qpole_k, :E_values, Float64[])) gamma=$(get(qpole_k, :gamma_values, Float64[]))")
        for scheme in (:phase_shift_current, :phase_shift_gbu_reference),
            convention in (:arg_propagator, :arg_inverse_propagator),
            anchor in (:none, :high_energy_zero),
            measure in (:single_charge_domega_over_pi, :legacy_domega_over_2pi)
            result = try
                phase_shift_meson_number_density(
                    pi_channel, qp, tp;
                    μ=μ_pi, degeneracy=1, scheme=scheme, qmax=4.0, q_nodes=4,
                    omega_min=0.4, omega_max=3.0, omega_nodes=8, eta=1e-3,
                    phase_convention=convention, phase_anchor=anchor,
                    omega_measure=measure, density_policy=:strict_normal_domain,
                )
            catch err
                (status=:error, density=NaN, message=sprint(showerror, err))
            end
            result_k = try
                phase_shift_meson_number_density(
                    k_channel, qp, tp;
                    μ=μ_K, degeneracy=1, scheme=scheme, qmax=4.0, q_nodes=4,
                    omega_min=0.4, omega_max=3.0, omega_nodes=8, eta=1e-3,
                    phase_convention=convention, phase_anchor=anchor,
                    omega_measure=measure, density_policy=:strict_normal_domain,
                )
            catch err
                (status=:error, density=NaN, message=sprint(showerror, err))
            end
            println(io, "phase scheme=$(scheme) convention=$(convention) anchor=$(anchor) measure=$(measure) " *
                "pi_status=$(result.status) K_status=$(result_k.status) " *
                "n_pi=$(_f(result.density)) n_K=$(_f(result_k.density)) " *
                "ratio=$(_f(result_k.density/result.density))")
        end
    end
end

function _phase_profiles(point, qp, tp, io)
    coeffs = _build_k_coeffs(qp)
    omega_grid, _ = gauleg(0.4, 3.0, 16)
    println(io, "[phase_profiles]")
    for meson in (:pi_plus, :pi_minus, :K_plus, :K_minus)
        for q in (0.0, 0.35)
            raw_prop = Float64[_propagator_phase(
                meson, ω, q, qp, tp, coeffs;
                eta=1e-3, real_axis_mode=:finite_eta,
                phase_convention=:arg_propagator,
            ) for ω in omega_grid]
            raw_inv = Float64[_propagator_phase(
                meson, ω, q, qp, tp, coeffs;
                eta=1e-3, real_axis_mode=:finite_eta,
                phase_convention=:arg_inverse_propagator,
            ) for ω in omega_grid]
            anchored = anchor_phase_high_energy(omega_grid, raw_prop; target=0.0)
            gbu = Float64[_phase_shift_weighted_phase(δ, :phase_shift_gbu_reference) for δ in anchored.anchored_phase]
            bose = Float64[bose_distribution(ω, 0.0, tp.T) * (1.0 + bose_distribution(ω, 0.0, tp.T)) for ω in omega_grid]
            println(io, "meson=$(meson) q=$(_f(q))")
            println(io, "omega=" * join(_f.(omega_grid), ","))
            println(io, "raw_argD=" * join(_f.(raw_prop), ","))
            println(io, "raw_minus_argInv=" * join(_f.(raw_inv), ","))
            println(io, "anchored_argD=" * join(_f.(anchored.anchored_phase), ","))
            println(io, "gbu_weighted=" * join(_f.(gbu), ","))
            println(io, "bose_g1g=" * join(_f.(bose), ","))
            println(io, "tail_before=$(_f(anchored.high_energy_phase_before_anchor)) shift=$(_f(anchored.applied_shift)) tail_span=$(_f(anchored.tail_span))")
            println(io, "ranges raw_argD=($(_f(minimum(raw_prop))),$(_f(maximum(raw_prop)))) raw_minus_argInv=($(_f(minimum(raw_inv))),$(_f(maximum(raw_inv)))) anchored=($(_f(minimum(anchored.anchored_phase))),$(_f(maximum(anchored.anchored_phase))))")
        end
    end
end

function _phase_shell_breakdown(point, qp, tp, io)
    coeffs = _build_k_coeffs(qp)
    q_grid, q_w = gauleg(0.0, 4.0, 4)
    omega_grid, omega_w = gauleg(0.4, 3.0, 8)
    println(io, "[phase_shell_breakdown]")
    for meson in (:K_plus, :K_minus)
        for convention in (:arg_propagator, :arg_inverse_propagator)
            μ = meson === :K_plus ? Float64(qp.μ.u - qp.μ.s) : Float64(qp.μ.s - qp.μ.u)
            total_none = 0.0
            total_anchor = 0.0
            println(io, "meson=$(meson) convention=$(convention) mu=$(_f(μ))")
            for iq in eachindex(q_grid, q_w)
                q = q_grid[iq]
                qw = q_w[iq]
                raw = Float64[_propagator_phase(meson, ω, q, qp, tp, coeffs;
                    eta=1e-3, real_axis_mode=:finite_eta, phase_convention=convention) for ω in omega_grid]
                profile = anchor_phase_high_energy(omega_grid, raw; target=0.0)
                g_none = 0.0
                g_anchor = 0.0
                for iω in eachindex(omega_grid, omega_w, raw, profile.anchored_phase)
                    ω = omega_grid[iω]
                    wω = omega_w[iω]
                    δraw = raw[iω]
                    δanch = profile.anchored_phase[iω]
                    g = bose_distribution(ω, μ, tp.T) * (1.0 + bose_distribution(ω, μ, tp.T))
                    g_none += wω * g * _phase_shift_weighted_phase(δraw, :phase_shift_gbu_reference)
                    g_anchor += wω * g * _phase_shift_weighted_phase(δanch, :phase_shift_gbu_reference)
                end
                shell_factor = q^2 / (2.0 * π^2) / tp.T
                shell_none = qw * shell_factor * g_none / π
                shell_anchor = qw * shell_factor * g_anchor / π
                total_none += shell_none
                total_anchor += shell_anchor
                println(io, "q[$iq]=$(_f(q)) raw_min=$(_f(minimum(raw))) raw_max=$(_f(maximum(raw))) " *
                    "anch_min=$(_f(minimum(profile.anchored_phase))) anch_max=$(_f(maximum(profile.anchored_phase))) " *
                    "tail_before=$(_f(profile.high_energy_phase_before_anchor)) shift=$(_f(profile.applied_shift)) " *
                    "raw_integral=$(_f(g_none)) anchored_integral=$(_f(g_anchor)) " *
                    "shell_none=$(_f(shell_none)) shell_anchor=$(_f(shell_anchor))")
            end
            println(io, "totals none=$(_f(total_none)) anchored=$(_f(total_anchor))")
        end
    end
end

function main()
    mkpath(dirname(OUTPUT_PATH))
    model, eq, point, qp, tp = _solve_state()
    open(OUTPUT_PATH, "w") do io
        println(io, "charged_rpa_bu_negative_density_diagnostic")
        println(io, "T_MeV=170 muB_MeV=240 Q_over_B=0.4 rhoS_target_fm3=0")
        println(io, "equilibrium residual=$(_f(eq.residual_norm)) converged=$(eq.converged)")
        println(io, "mu_u=$(_f(qp.μ.u)) mu_d=$(_f(qp.μ.d)) mu_s=$(_f(qp.μ.s))")
        println(io, "m_u=$(_f(qp.m.u)) m_d=$(_f(qp.m.d)) m_s=$(_f(qp.m.s))")
        println(io, "m_pi=$(_f(point.meson_results[:pi].mass)) gamma_pi=$(_f(point.meson_results[:pi].gamma))")
        println(io, "m_K=$(_f(point.meson_results[:K].mass)) gamma_K=$(_f(point.meson_results[:K].gamma))")
        println(io)
        _density_suite(point, qp, tp, io)
        println(io)
        _phase_shell_breakdown(point, qp, tp, io)
        println(io)
        _phase_profiles(point, qp, tp, io)
    end
    println("wrote ", OUTPUT_PATH)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
