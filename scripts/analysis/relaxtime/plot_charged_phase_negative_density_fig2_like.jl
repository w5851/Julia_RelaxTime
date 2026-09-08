"""
Analysis-only visualization of the strict charged phase/BU negative-density
diagnostic.

The script uses one finite-BQS quark-only background and the same ordered
charged-RPA inverse propagator as `audit_charged_phase_backend.jl`.  It keeps
the phase object explicit:

    delta(omega,q) = -arg(Delta_inverse^R(omega,q))

For every channel and sampled q it writes the principal (raw), continuous
unwrapped, and high-energy-anchored phases, together with both the current BU
weight and the Fig.2 generalized-BU weight
`F(delta) = delta - sin(2*delta)/2`.

The two figures are diagnostic only.  They do not change a production default,
the PNJL equilibrium solver, or any density baseline.
"""

const AUDIT_SCRIPT = joinpath(@__DIR__, "audit_charged_phase_backend.jl")
include(AUDIT_SCRIPT)

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using CSV
using JSON3
using SHA
using Plots
using Printf: @sprintf
using Main.Constants_PNJL: ħc_MeV_fm
using Main.RelaxTime.GaussLegendre: gauleg
using Main.RelaxTime.AFieldBuilder: build_A_triplet
using Main.RelaxTime.ChargedRPAKernel: charged_rpa_spec, charged_rpa_coupling,
                                       charged_rpa_inverse
using Main.RelaxTime.ChargedRPAProvider: charged_polarization,
                                          charged_pair_continuum_thresholds
using Main.RelaxTime.ChargedPhaseBackend: StrictChargedPhaseSpec,
                                          strict_phase_profile
using Main.RelaxTime.MesonInteractionKernel: build_full_kmt_interaction

const DEFAULT_OUTPUT_DIR = joinpath(
    PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "analysis",
    "charged_rpa_phase_backend", "negative_density_phase_fig2_like",
)

const CHANNEL_LABELS = Dict(
    :pi_plus => "pi+",
    :pi_minus => "pi-",
    :K_plus => "K+",
    :K_minus => "K-",
)

const CHANNEL_COLORS = Dict(
    :pi_plus => :royalblue,
    :pi_minus => :darkorange,
    :K_plus => :forestgreen,
    :K_minus => :firebrick,
)

@inline _plot_float(name::AbstractString, default::Real) =
    parse(Float64, get(ENV, "CHARGED_PHASE_PLOT_" * name, string(default)))
@inline _plot_int(name::AbstractString, default::Integer) =
    parse(Int, get(ENV, "CHARGED_PHASE_PLOT_" * name, string(default)))

@inline function _legacy_env_int(name::AbstractString, default::Integer)
    return parse(Int, get(ENV, name, string(default)))
end

@inline function _legacy_env_float(name::AbstractString, default::Real)
    return parse(Float64, get(ENV, name, string(default)))
end

function _parse_float_list(raw::AbstractString)
    values = Float64[]
    for piece in split(raw, ',')
        text = strip(piece)
        isempty(text) && continue
        push!(values, parse(Float64, text))
    end
    isempty(values) && throw(ArgumentError("a q/omega list cannot be empty"))
    return values
end

function _plot_q_values(qmax::Float64)
    raw = strip(get(ENV, "CHARGED_PHASE_PLOT_Q_VALUES", ""))
    values = isempty(raw) ? qmax .* [0.0, 0.25, 0.5, 0.75, 1.0] : _parse_float_list(raw)
    all(isfinite, values) && all(value -> 0.0 <= value <= qmax, values) ||
        throw(ArgumentError("CHARGED_PHASE_PLOT_Q_VALUES must lie in [0,qmax]"))
    return unique(sort(values))
end

function _plot_omega_grid(mu_values, T_fm)
    # Use one common grid for all four channels.  The lower edge is above the
    # largest charged chemical potential, so the normal-phase Bose factor is
    # defined for every plotted channel.
    margin = _plot_float("BOSE_MARGIN_INV_FM", 2.0e-3)
    omega_min = max(_plot_float("OMEGA_MIN_INV_FM", 0.05), maximum(mu_values) + margin)
    omega_max = _plot_float("OMEGA_MAX_INV_FM", 8.0)
    points = _plot_int("OMEGA_POINTS", 180)
    omega_max > omega_min || throw(ArgumentError("omega max must exceed omega min"))
    points >= 16 || throw(ArgumentError("omega points must be at least 16"))
    return collect(range(omega_min, omega_max; length=points))
end

@inline function _phase_derivative(x::Vector{Float64}, y::Vector{Float64})
    result = similar(y)
    result[1] = (y[2] - y[1]) / (x[2] - x[1])
    result[end] = (y[end] - y[end - 1]) / (x[end] - x[end - 1])
    @inbounds for i in 2:(length(y) - 1)
        result[i] = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
    end
    return result
end

@inline function _trapz(x::Vector{Float64}, y::Vector{Float64})
    total = 0.0
    @inbounds for i in 1:(length(x) - 1)
        total += 0.5 * (x[i + 1] - x[i]) * (y[i + 1] + y[i])
    end
    return total
end

@inline _gbu_phase(delta::Real) = Float64(delta) - 0.5 * sin(2.0 * Float64(delta))
@inline _gbu_phase_derivative(delta::Real) = 2.0 * sin(Float64(delta))^2

# Fig.2 uses a display-only representative in [0,pi].  Never feed this fold
# back into the BU derivative: the integration branch remains `anchored`.
@inline function _fold_0_pi(delta::Real)
    value = Float64(delta)
    return π - abs(mod(value, 2.0π) - π)
end

function _bose_weights(omega::Vector{Float64}, mu::Float64, T::Float64)
    bose = fill(NaN, length(omega))
    bose_g1g = fill(NaN, length(omega))
    support = fill("unsafe_omega_le_mu", length(omega))
    @inbounds for i in eachindex(omega)
        x = (omega[i] - mu) / T
        if x > 0.0
            g = x > 700.0 ? 0.0 : inv(expm1(x))
            bose[i] = g
            bose_g1g[i] = g * (1.0 + g)
            support[i] = "safe_normal_domain"
        end
    end
    return bose, bose_g1g, support
end

function _profile_spec()
    return StrictChargedPhaseSpec(
        phase_object=:inverse_propagator,
        phase_sign=-1,
        target=_plot_float("ANCHOR_TARGET", 0.0),
        branch_tol=_plot_float("BRANCH_TOL", 0.0),
        tail_points=_plot_int("TAIL_POINTS", 4),
        tail_tolerance=_plot_float("TAIL_TOLERANCE", 0.2π),
    )
end

function _variants()
    variants = [(name="pv_cut", prescription=:ordered_pv_cut, eta=0.0)]
    include_retarded = lowercase(get(ENV, "CHARGED_PHASE_PLOT_INCLUDE_RETARDED", "true")) in
        ("1", "true", "yes")
    if include_retarded
        eta = _plot_float("RETARDED_ETA_INV_FM", 5.0e-3)
        eta > 0.0 || throw(ArgumentError("retarded eta must be positive"))
        push!(variants, (name="retarded_eta_$(replace(@sprintf("%.4g", eta), "." => "p"))",
                         prescription=:ordered_retarded, eta=eta))
    end
    return variants
end

function _background_for_plot()
    background = _solve_background()
    state = Models.meanfield_state(background.result.x_state)
    masses = (
        u=Float64(background.result.masses[1]),
        d=Float64(background.result.masses[2]),
        s=Float64(background.result.masses[3]),
    )
    chemical_potentials = (
        u=Float64(background.result.mu_vec[1]),
        d=Float64(background.result.mu_vec[2]),
        s=Float64(background.result.mu_vec[3]),
    )
    thermo = (
        T=background.T_fm,
        Φ=Float64(state.Phi),
        Φbar=Float64(state.PhiBar),
        ξ=0.0,
    )
    A_values = build_A_triplet(
        (m=masses, μ=chemical_potentials),
        thermo;
        p_nodes=_plot_int("A_NODES", _legacy_env_int("CHARGED_PHASE_A_NODES", 64)),
        p_max=_plot_float("A_PMAX", _legacy_env_float("CHARGED_PHASE_A_PMAX", 16.0)),
        use_aniso=false,
    )
    kernel = build_full_kmt_interaction(
        state.phi;
        G=background.model.params.G_fm2,
        K=background.model.params.K_fm5,
    )
    return background, masses, chemical_potentials, thermo, A_values, kernel
end

function _channel_mass_and_mu(channel::Symbol, background, chemical_potentials)
    if channel === :pi_plus || channel === :pi_minus
        mass = Float64(background.point.meson_results[:pi].mass)
        mu = channel === :pi_plus ? chemical_potentials.u - chemical_potentials.d :
            chemical_potentials.d - chemical_potentials.u
    elseif channel === :K_plus || channel === :K_minus
        mass = Float64(background.point.meson_results[:K].mass)
        mu = channel === :K_plus ? chemical_potentials.u - chemical_potentials.s :
            chemical_potentials.s - chemical_potentials.u
    else
        throw(ArgumentError("unsupported charged channel $(channel)"))
    end
    return mass, Float64(mu)
end

function _inverse_builder(channel, masses, chemical_potentials, thermo, A_values, kernel, variant)
    spec = charged_rpa_spec(channel)
    coupling = charged_rpa_coupling(kernel, spec)
    polarization(omega, q) = charged_polarization(
        spec,
        omega,
        q,
        masses,
        chemical_potentials,
        thermo,
        A_values;
        prescription=variant.prescription,
        eta_inv_fm=variant.eta,
        energy_nodes=_plot_int("POLARIZATION_NODES", 24),
    ).value
    inverse(omega, q) = charged_rpa_inverse(spec, coupling, polarization(omega, q))
    return spec, Float64(coupling), inverse
end

function _evaluate_profile(channel, q, omega, background, masses, chemical_potentials,
                           thermo, A_values, kernel, variant, phase_spec)
    spec, coupling, inverse = _inverse_builder(
        channel, masses, chemical_potentials, thermo, A_values, kernel, variant,
    )
    inverse_values = ComplexF64[]
    error_message = ""
    for value in omega
        try
            push!(inverse_values, ComplexF64(inverse(value, q)))
        catch err
            error_message = sprint(showerror, err)
            break
        end
    end
    if length(inverse_values) != length(omega)
        return (
            ok=false,
            error_message=error_message,
            profile=nothing,
            inverse_values=ComplexF64[],
            coupling=coupling,
            spec=spec,
        )
    end
    profile = try
        strict_phase_profile(omega, inverse_values; spec=phase_spec)
    catch err
        return (
            ok=false,
            error_message=sprint(showerror, err),
            profile=nothing,
            inverse_values=inverse_values,
            coupling=coupling,
            spec=spec,
        )
    end
    return (
        ok=true,
        error_message="",
        profile=profile,
        inverse_values=inverse_values,
        coupling=coupling,
        spec=spec,
    )
end

function _shell_summary(profile, omega, bose, bose_g1g, q, T_fm)
    anchored = Float64.(profile.anchored_phase)
    unwrapped = Float64.(profile.unwrapped_phase)
    raw = Float64.(profile.raw_phase)
    d_anchored = _phase_derivative(omega, anchored)
    d_unwrapped = _phase_derivative(omega, unwrapped)
    d_raw = _phase_derivative(omega, raw)
    # Difference F itself: the smooth chain rule loses unresolved pi jumps.
    gbu_derivative = _phase_derivative(omega, _gbu_phase.(anchored))
    gbu_derivative_raw = _phase_derivative(omega, _gbu_phase.(raw))
    current_integrand = bose .* d_anchored
    current_integrand_raw = bose .* d_raw
    current_integrand_unwrapped = bose .* d_unwrapped
    gbu_integrand = bose .* gbu_derivative
    gbu_integrand_raw = bose .* gbu_derivative_raw
    current_g1g_integrand = bose_g1g .* anchored
    q_prefactor = q^2 / (2.0 * π^2)
    measure_factor = inv(π)
    current_shell = q_prefactor * measure_factor * _trapz(omega, current_integrand)
    current_shell_raw = q_prefactor * measure_factor * _trapz(omega, current_integrand_raw)
    current_shell_unwrapped = q_prefactor * measure_factor * _trapz(omega, current_integrand_unwrapped)
    gbu_shell = q_prefactor * measure_factor * _trapz(omega, gbu_integrand)
    gbu_shell_raw = q_prefactor * measure_factor * _trapz(omega, gbu_integrand_raw)
    return (
        raw=raw,
        unwrapped=unwrapped,
        anchored=anchored,
        d_raw=d_raw,
        d_unwrapped=d_unwrapped,
        d_anchored=d_anchored,
        gbu=Float64[_gbu_phase(value) for value in anchored],
        phase_display=Float64[_fold_0_pi(value) for value in anchored],
        gbu_display=Float64[_gbu_phase(_fold_0_pi(value)) for value in anchored],
        current_integrand=current_integrand,
        current_integrand_raw=current_integrand_raw,
        current_integrand_unwrapped=current_integrand_unwrapped,
        gbu_integrand=gbu_integrand,
        gbu_integrand_raw=gbu_integrand_raw,
        current_g1g_integrand=current_g1g_integrand,
        current_shell_fm3_per_dq=current_shell,
        current_shell_raw_fm3_per_dq=current_shell_raw,
        current_shell_unwrapped_fm3_per_dq=current_shell_unwrapped,
        gbu_shell_fm3_per_dq=gbu_shell,
        gbu_shell_raw_fm3_per_dq=gbu_shell_raw,
        min_current_integrand=minimum(current_integrand),
        max_current_integrand=maximum(current_integrand),
        min_current_integrand_raw=minimum(current_integrand_raw),
        max_current_integrand_raw=maximum(current_integrand_raw),
        negative_current_fraction=count(<(0.0), current_integrand) / length(current_integrand),
        negative_current_raw_fraction=count(<(0.0), current_integrand_raw) / length(current_integrand_raw),
        negative_gbu_fraction=count(<(0.0), gbu_integrand) / length(gbu_integrand),
        negative_gbu_raw_fraction=count(<(0.0), gbu_integrand_raw) / length(gbu_integrand_raw),
    )
end

function _failure_shell(omega)
    nan_values = fill(NaN, length(omega))
    return (
        raw=nan_values,
        unwrapped=nan_values,
        anchored=nan_values,
        d_raw=nan_values,
        d_unwrapped=nan_values,
        d_anchored=nan_values,
        gbu=nan_values,
        phase_display=nan_values,
        gbu_display=nan_values,
        current_integrand=nan_values,
        current_integrand_raw=nan_values,
        current_integrand_unwrapped=nan_values,
        gbu_integrand=nan_values,
        gbu_integrand_raw=nan_values,
        current_g1g_integrand=nan_values,
        current_shell_fm3_per_dq=NaN,
        current_shell_raw_fm3_per_dq=NaN,
        current_shell_unwrapped_fm3_per_dq=NaN,
        gbu_shell_fm3_per_dq=NaN,
        gbu_shell_raw_fm3_per_dq=NaN,
        min_current_integrand=NaN,
        max_current_integrand=NaN,
        min_current_integrand_raw=NaN,
        max_current_integrand_raw=NaN,
        negative_current_fraction=NaN,
        negative_current_raw_fraction=NaN,
        negative_gbu_fraction=NaN,
        negative_gbu_raw_fraction=NaN,
    )
end

function _phase_rows!(rows, shell_rows, channel, q, omega, profile_result, shell,
                      bose, bose_g1g, support, background, masses, chemical_potentials,
                      variant, threshold, mass, mu, thermo)
    profile = profile_result.profile
    inverse_values = profile_result.inverse_values
    if profile_result.ok
        phase_status = "ok"
        tail_stable = Bool(profile.tail_stable)
        anchor_shift = Float64(profile.applied_shift)
        tail_span = Float64(profile.tail_span)
        high_energy_before = Float64(profile.high_energy_phase_before_anchor)
        high_energy_after = Float64(profile.high_energy_phase_after_anchor)
    else
        phase_status = "evaluation_error:$(replace(profile_result.error_message, ',' => ';'))"
        tail_stable = false
        anchor_shift = NaN
        tail_span = NaN
        high_energy_before = NaN
        high_energy_after = NaN
    end
    threshold_mev = threshold * ħc_MeV_fm
    for i in eachindex(omega)
        inv_value = profile_result.ok ? inverse_values[i] : ComplexF64(NaN, NaN)
        push!(rows, (
            channel=String(channel),
            channel_label=CHANNEL_LABELS[channel],
            variant=variant.name,
            prescription=String(variant.prescription),
            phase_object="inverse_propagator",
            phase_sign=-1,
            anchor="high_energy_zero",
            T_MeV=background.T_MeV,
            muB_MeV=background.muB_MeV,
            q_inv_fm=Float64(q),
            omega_inv_fm=Float64(omega[i]),
            omega_MeV=Float64(omega[i] * ħc_MeV_fm),
            threshold_inv_fm=Float64(threshold),
            threshold_MeV=Float64(threshold_mev),
            meson_mass_inv_fm=Float64(mass),
            meson_mu_inv_fm=Float64(mu),
            meson_mu_MeV=Float64(mu * ħc_MeV_fm),
            mu_u_inv_fm=Float64(chemical_potentials.u),
            mu_d_inv_fm=Float64(chemical_potentials.d),
            mu_s_inv_fm=Float64(chemical_potentials.s),
            m_u_inv_fm=Float64(masses.u),
            m_d_inv_fm=Float64(masses.d),
            m_s_inv_fm=Float64(masses.s),
            coupling_fm2=Float64(profile_result.coupling),
            inverse_real=real(inv_value),
            inverse_imag=imag(inv_value),
            s_matrix_real=profile_result.ok ? real(profile.s_matrix_values[i]) : NaN,
            s_matrix_imag=profile_result.ok ? imag(profile.s_matrix_values[i]) : NaN,
            raw_phase=Float64(shell.raw[i]),
            unwrapped_phase=Float64(shell.unwrapped[i]),
            anchored_phase=Float64(shell.anchored[i]),
            gbu_phase=Float64(shell.gbu[i]),
            phase_display_fold_0_pi=Float64(shell.phase_display[i]),
            gbu_display_fold_0_pi=Float64(shell.gbu_display[i]),
            phase_derivative_raw=Float64(shell.d_raw[i]),
            phase_derivative_unwrapped=Float64(shell.d_unwrapped[i]),
            phase_derivative_anchored=Float64(shell.d_anchored[i]),
            bose_weight=Float64(bose[i]),
            bose_g1g_weight=Float64(bose_g1g[i]),
            current_integrand=Float64(shell.current_integrand[i]),
            current_integrand_raw=Float64(shell.current_integrand_raw[i]),
            current_integrand_unwrapped=Float64(shell.current_integrand_unwrapped[i]),
            gbu_integrand=Float64(shell.gbu_integrand[i]),
            gbu_integrand_raw=Float64(shell.gbu_integrand_raw[i]),
            current_g1g_integrand=Float64(shell.current_g1g_integrand[i]),
            omega_measure="single_charge_domega_over_pi",
            omega_measure_factor=inv(π),
            bose_support_status=support[i],
            q_shell_current_fm3_per_dq=Float64(shell.current_shell_fm3_per_dq),
            q_shell_current_raw_fm3_per_dq=Float64(shell.current_shell_raw_fm3_per_dq),
            q_shell_current_unwrapped_fm3_per_dq=Float64(shell.current_shell_unwrapped_fm3_per_dq),
            q_shell_gbu_fm3_per_dq=Float64(shell.gbu_shell_fm3_per_dq),
            q_shell_gbu_raw_fm3_per_dq=Float64(shell.gbu_shell_raw_fm3_per_dq),
            q_shell_current_negative=Bool(isfinite(shell.current_shell_fm3_per_dq) && shell.current_shell_fm3_per_dq < 0.0),
            q_shell_gbu_negative=Bool(isfinite(shell.gbu_shell_fm3_per_dq) && shell.gbu_shell_fm3_per_dq < 0.0),
            phase_status=phase_status,
            tail_stable=tail_stable,
            anchor_shift=anchor_shift,
            tail_span=tail_span,
            high_energy_phase_before_anchor=high_energy_before,
            high_energy_phase_after_anchor=high_energy_after,
        ))
    end
    push!(shell_rows, (
        channel=String(channel),
        channel_label=CHANNEL_LABELS[channel],
        variant=variant.name,
        prescription=String(variant.prescription),
        q_inv_fm=Float64(q),
        q_shell_current_fm3_per_dq=Float64(shell.current_shell_fm3_per_dq),
        q_shell_current_raw_fm3_per_dq=Float64(shell.current_shell_raw_fm3_per_dq),
        q_shell_current_unwrapped_fm3_per_dq=Float64(shell.current_shell_unwrapped_fm3_per_dq),
        q_shell_gbu_fm3_per_dq=Float64(shell.gbu_shell_fm3_per_dq),
        q_shell_gbu_raw_fm3_per_dq=Float64(shell.gbu_shell_raw_fm3_per_dq),
        min_current_integrand=Float64(shell.min_current_integrand),
        max_current_integrand=Float64(shell.max_current_integrand),
        min_current_integrand_raw=Float64(shell.min_current_integrand_raw),
        max_current_integrand_raw=Float64(shell.max_current_integrand_raw),
        negative_current_fraction=Float64(shell.negative_current_fraction),
        negative_current_raw_fraction=Float64(shell.negative_current_raw_fraction),
        negative_gbu_fraction=Float64(shell.negative_gbu_fraction),
        negative_gbu_raw_fraction=Float64(shell.negative_gbu_raw_fraction),
        phase_status=phase_status,
        tail_stable=tail_stable,
        threshold_MeV=Float64(threshold_mev),
        meson_mu_MeV=Float64(mu * ħc_MeV_fm),
        T_MeV=background.T_MeV,
        muB_MeV=background.muB_MeV,
        gap_residual_norm=Float64(background.result.residual_norm),
        omega_min_MeV=Float64(first(omega) * ħc_MeV_fm),
        omega_max_MeV=Float64(last(omega) * ħc_MeV_fm),
        omega_points=length(omega),
        polarization_nodes=_plot_int("POLARIZATION_NODES", 24),
    ))
end

function _select_rows(rows, channel::Symbol, variant::AbstractString, q::Real)
    selected = [row for row in rows if row.channel == String(channel) &&
        row.variant == variant && isapprox(row.q_inv_fm, Float64(q); atol=1.0e-10)]
    sort!(selected; by=row -> row.omega_inv_fm)
    return selected
end

function _finite_xy(rows, yfield::Symbol)
    filtered = [row for row in rows if isfinite(Float64(getproperty(row, :omega_MeV))) &&
        isfinite(Float64(getproperty(row, yfield)))]
    return ([Float64(row.omega_MeV) / 1000.0 for row in filtered],
            [Float64(getproperty(row, yfield)) for row in filtered])
end

function _render_fig2_like(rows, q_values, variants, output_path, background)
    pv_name = first(variants).name
    q0 = first(q_values)
    qmid = q_values[cld(length(q_values), 2)]
    x_plot_max = _plot_float("XMAX_GEV", 2.2)
    p_left = plot(; xlabel="omega [GeV]", ylabel="delta(omega,q)",
        title="Finite-BQS fold_0_pi display, q=$(q0)", legend=:outertopright, grid=true,
        xlims=(0.0, x_plot_max), size=(520, 430), linewidth=2.0)
    for channel in CHARGED_MODES
        selected = _select_rows(rows, channel, pv_name, q0)
        x, y = _finite_xy(selected, :phase_display_fold_0_pi)
        plot!(p_left, x, y; label=CHANNEL_LABELS[channel], color=CHANNEL_COLORS[channel])
        if !isempty(selected)
            threshold = Float64(first(selected).threshold_MeV) / 1000.0
            vline!(p_left, [threshold]; label="", color=CHANNEL_COLORS[channel],
                linestyle=:dash, alpha=0.45, linewidth=0.9)
        end
    end
    plot!(p_left, [0.0, x_plot_max], [0.0, 0.0]; label="", color=:black, linestyle=:dot, linewidth=0.8)

    p_center = plot(; xlabel="omega [GeV]", ylabel="delta(omega,q)",
        title="K+ phase branch comparison", legend=:outertopright, grid=true,
        xlims=(0.0, x_plot_max), size=(520, 430), linewidth=2.0)
    center_specs = [(pv_name, q0, "PV folded display q=$(q0)", :forestgreen, :solid),
                    (pv_name, q0, "PV anchored branch q=$(q0)", :forestgreen, :dash),
                    (pv_name, qmid, "PV folded display q=$(round(qmid; digits=3))", :darkgreen, :dot)]
    if length(variants) > 1
        ret_name = variants[2].name
        push!(center_specs, (ret_name, q0, "retarded folded display q=$(q0)", :crimson, :solid))
        push!(center_specs, (ret_name, q0, "retarded anchored branch q=$(q0)", :crimson, :dash))
    end
    for (variant_name, q, label, color, style) in center_specs
        selected = _select_rows(rows, :K_plus, variant_name, q)
        field = occursin("anchored branch", label) ? :anchored_phase : :phase_display_fold_0_pi
        x, y = _finite_xy(selected, field)
        plot!(p_center, x, y; label=label, color=color, linestyle=style)
    end
    kplus_q0 = _select_rows(rows, :K_plus, pv_name, q0)
    if !isempty(kplus_q0)
        vline!(p_center, [Float64(first(kplus_q0).threshold_MeV) / 1000.0]; label="K+ threshold",
            color=:black, linestyle=:dash, alpha=0.55, linewidth=0.9)
    end
    plot!(p_center, [0.0, x_plot_max], [0.0, 0.0]; label="", color=:black, linestyle=:dot, linewidth=0.8)

    p_right = plot(; xlabel="omega [GeV]", ylabel="delta - sin(2delta)/2",
        title="Fig.2 generalized-BU weight", legend=:outertopright, grid=true,
        xlims=(0.0, x_plot_max), size=(520, 430), linewidth=2.0)
    for (variant_name, q, label, color, style) in center_specs
        selected = _select_rows(rows, :K_plus, variant_name, q)
        field = occursin("anchored branch", label) ? :gbu_phase : :gbu_display_fold_0_pi
        x, y = _finite_xy(selected, field)
        plot!(p_right, x, y; label=label, color=color, linestyle=style)
    end
    if !isempty(kplus_q0)
        vline!(p_right, [Float64(first(kplus_q0).threshold_MeV) / 1000.0]; label="K+ threshold",
            color=:black, linestyle=:dash, alpha=0.55, linewidth=0.9)
    end
    plot!(p_right, [0.0, x_plot_max], [0.0, 0.0]; label="", color=:black, linestyle=:dot, linewidth=0.8)

    fig = plot(p_left, p_center, p_right; layout=(1, 3),
        size=(1640, 480), margin=5Plots.mm,
        plot_title="Charged-RPA phase diagnostic: T=$(background.T_MeV) MeV, muB=$(background.muB_MeV) MeV")
    savefig(fig, output_path)
end

function _render_attribution(rows, shell_rows, q_values, variants, output_path, background)
    pv_name = first(variants).name
    q0 = first(q_values)
    x_plot_max = _plot_float("XMAX_GEV", 2.2)
    p_integrand = plot(; xlabel="omega [GeV]", ylabel="g_B * ddelta/domega",
        title="PV current BU integrand at q=$(q0)", legend=:outertopright, grid=true,
        xlims=(0.0, x_plot_max), size=(720, 430), linewidth=2.0)
    for channel in CHARGED_MODES
        selected = _select_rows(rows, channel, pv_name, q0)
        x, y = _finite_xy(selected, :current_integrand)
        plot!(p_integrand, x, y; label=CHANNEL_LABELS[channel], color=CHANNEL_COLORS[channel],
            fillrange=0.0, fillalpha=0.08)
        x_raw, y_raw = _finite_xy(selected, :current_integrand_raw)
        plot!(p_integrand, x_raw, y_raw; label="", color=CHANNEL_COLORS[channel], linestyle=:dash,
            linewidth=1.0)
    end
    plot!(p_integrand, [0.0, x_plot_max], [0.0, 0.0]; label="zero", color=:black, linestyle=:dot, linewidth=0.9)

    p_shell = plot(; xlabel="q [fm^-1]", ylabel="shell density / dq [fm^-3 per fm^-1]",
        title="PV q-shell sign (solid=anchored, dash=raw, dot=GBU)", legend=:outertopright, grid=true,
        size=(720, 430), linewidth=2.0)
    for channel in CHARGED_MODES
        selected = [row for row in shell_rows if row.channel == String(channel) && row.variant == pv_name]
        sort!(selected; by=row -> row.q_inv_fm)
        q = [Float64(row.q_inv_fm) for row in selected]
        current = [Float64(row.q_shell_current_fm3_per_dq) for row in selected]
        current_raw = [Float64(row.q_shell_current_raw_fm3_per_dq) for row in selected]
        gbu = [Float64(row.q_shell_gbu_fm3_per_dq) for row in selected]
        color = CHANNEL_COLORS[channel]
        plot!(p_shell, q, current; label=CHANNEL_LABELS[channel], color=color, linestyle=:solid)
        plot!(p_shell, q, current_raw; label="", color=color, linestyle=:dash, linewidth=1.0)
        plot!(p_shell, q, gbu; label="", color=color, linestyle=:dot, linewidth=1.0)
    end
    plot!(p_shell, [minimum(q_values), maximum(q_values)], [0.0, 0.0]; label="zero", color=:black, linestyle=:dot, linewidth=0.9)

    fig = plot(p_integrand, p_shell; layout=(1, 2), size=(1480, 500), margin=5Plots.mm,
        plot_title="Negative-density attribution: PV phase derivative and q shells")
    savefig(fig, output_path)
end

function _write_manifest(path, script_path, output_dir, background, masses, chemical_potentials,
                         thermo, q_values, omega, variants, phase_spec, rows, shell_rows)
    head = try
        readchomp(`git -C $(PROJECT_ROOT) rev-parse HEAD`)
    catch
        "unknown"
    end
    manifest = Dict(
        "script" => relpath(script_path, PROJECT_ROOT),
        "git_head" => head,
        "schema" => "charged_phase_plot_v2",
        "source_sha256" => bytes2hex(sha256(read(script_path))),
        "worktree_dirty" => !isempty(readchomp(`git -C $(PROJECT_ROOT) status --porcelain --untracked-files=no`)),
        "threshold_coordinate" => "external_k0",
        "threshold_formula" => "hypot(q,m1+m2)-(mu1-mu2)",
        "shell_units" => "fm^-2 = fm^-3 per fm^-1",
        "shell_measure" => "g_B dF/pi; no extra 1/T",
        "bose_coordinate_status" => "legacy g(k0-mu_M) diagnostic; KMS coordinate review required",
        "root_count_status" => "unwrap is not independent bound-state counting",
        "route" => "analysis_only_strict_charged_phase_negative_density",
        "phase_object" => "inverse_propagator",
        "phase_definition" => "delta = -arg(Delta_inverse^R)",
        "phase_weight" => "delta - sin(2 delta)/2",
        "phase_display" => "fold_0_pi = pi - abs(mod(delta,2pi)-pi); display only",
        "background" => Dict(
            "T_MeV" => background.T_MeV,
            "muB_MeV" => background.muB_MeV,
            "Q_over_B" => 0.4,
            "rhoS_target" => 0.0,
            "gap_residual_norm" => background.result.residual_norm,
        ),
        "masses_inv_fm" => Dict("u" => masses.u, "d" => masses.d, "s" => masses.s),
        "chemical_potentials_inv_fm" => Dict("u" => chemical_potentials.u,
                                             "d" => chemical_potentials.d,
                                             "s" => chemical_potentials.s),
        "thermo" => Dict("T_inv_fm" => thermo.T, "Phi" => thermo.Φ, "PhiBar" => thermo.Φbar),
        "q_values_inv_fm" => q_values,
        "omega_grid_inv_fm" => Dict("min" => first(omega), "max" => last(omega), "points" => length(omega)),
        "variants" => [Dict("name" => v.name, "prescription" => String(v.prescription), "eta_inv_fm" => v.eta) for v in variants],
        "phase_spec" => Dict("target" => phase_spec.target, "tail_points" => phase_spec.tail_points,
                             "tail_tolerance" => phase_spec.tail_tolerance, "branch_tol" => phase_spec.branch_tol),
        "rows" => length(rows),
        "A_nodes" => _plot_int("A_NODES", _legacy_env_int("CHARGED_PHASE_A_NODES", 64)),
        "A_pmax_inv_fm" => _plot_float("A_PMAX", _legacy_env_float("CHARGED_PHASE_A_PMAX", 16.0)),
        "requested_polarization_nodes" => _plot_int("POLARIZATION_NODES", 24),
        "pv_quadrature" => "OneLoopIntegrals fixed hybrid 16/32; requested polarization nodes apply to finite eta only",
        "shell_rows" => length(shell_rows),
        "outputs" => Dict(
            "detail_csv" => "charged_phase_profile_detail.csv",
            "q_shell_csv" => "charged_phase_q_shell_summary.csv",
            "fig2_like" => "charged_phase_fig2_like.png",
            "attribution" => "charged_phase_negative_density_attribution.png",
        ),
        "status" => "diagnostic_only_not_production",
    )
    open(path, "w") do io
        JSON3.write(io, manifest)
    end
end

function main()
    output_dir = abspath(get(ENV, "CHARGED_PHASE_PLOT_OUTPUT_DIR", DEFAULT_OUTPUT_DIR))
    isdir(output_dir) && !isempty(readdir(output_dir)) &&
        error("refusing to overwrite existing diagnostic output directory: $(output_dir)")
    mkpath(output_dir)
    background, masses, chemical_potentials, thermo, A_values, kernel = _background_for_plot()
    channels = CHARGED_MODES
    channel_mu = [_channel_mass_and_mu(channel, background, chemical_potentials)[2] for channel in channels]
    omega = _plot_omega_grid(channel_mu, thermo.T)
    qmax = _plot_float("QMAX_INV_FM", 1.0)
    qmax > 0.0 || throw(ArgumentError("plot qmax must be positive"))
    q_values = _plot_q_values(qmax)
    variants = _variants()
    phase_spec = _profile_spec()
    rows = NamedTuple[]
    shell_rows = NamedTuple[]
    for variant in variants
        for channel in channels
            mass, mu = _channel_mass_and_mu(channel, background, chemical_potentials)
            spec = charged_rpa_spec(channel)
            m1 = getproperty(masses, spec.pair[1])
            m2 = getproperty(masses, spec.pair[2])
            mu1 = getproperty(chemical_potentials, spec.pair[1])
            mu2 = getproperty(chemical_potentials, spec.pair[2])
            for q in q_values
                threshold = charged_pair_continuum_thresholds(q, m1, m2, mu1, mu2).k0_threshold_inv_fm
                profile_result = _evaluate_profile(
                    channel, q, omega, background, masses, chemical_potentials,
                    thermo, A_values, kernel, variant, phase_spec,
                )
                bose, bose_g1g, support = _bose_weights(omega, mu, thermo.T)
                shell = profile_result.ok ? _shell_summary(profile_result.profile, omega, bose, bose_g1g, q, thermo.T) : _failure_shell(omega)
                _phase_rows!(
                    rows, shell_rows, channel, q, omega, profile_result, shell,
                    bose, bose_g1g, support, background, masses, chemical_potentials,
                    variant, threshold, mass, mu, thermo,
                )
                println(@sprintf("[charged-phase-plot] %-8s %-24s q=%6.3f status=%s shell_current=% .6e shell_gbu=% .6e",
                    CHANNEL_LABELS[channel], variant.name, q,
                    profile_result.ok ? "ok" : "failed",
                    shell.current_shell_fm3_per_dq, shell.gbu_shell_fm3_per_dq))
            end
        end
    end

    detail_path = joinpath(output_dir, "charged_phase_profile_detail.csv")
    shell_path = joinpath(output_dir, "charged_phase_q_shell_summary.csv")
    fig2_path = joinpath(output_dir, "charged_phase_fig2_like.png")
    attribution_path = joinpath(output_dir, "charged_phase_negative_density_attribution.png")
    manifest_path = joinpath(output_dir, "plot_manifest.json")
    CSV.write(detail_path, rows)
    CSV.write(shell_path, shell_rows)
    _render_fig2_like(rows, q_values, variants, fig2_path, background)
    _render_attribution(rows, shell_rows, q_values, variants, attribution_path, background)
    _write_manifest(
        manifest_path, @__FILE__, output_dir, background, masses, chemical_potentials,
        thermo, q_values, omega, variants, phase_spec, rows, shell_rows,
    )
    println("[charged-phase-plot] detail: $(detail_path)")
    println("[charged-phase-plot] q shells: $(shell_path)")
    println("[charged-phase-plot] Fig.2-like: $(fig2_path)")
    println("[charged-phase-plot] attribution: $(attribution_path)")
    return (rows=rows, shell_rows=shell_rows, detail_path=detail_path,
            shell_path=shell_path, fig2_path=fig2_path, attribution_path=attribution_path)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
