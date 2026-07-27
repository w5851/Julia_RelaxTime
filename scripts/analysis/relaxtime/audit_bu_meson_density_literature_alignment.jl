"""
BU 介子数密度文献对齐审查脚本。

目标：
- 固定最小 fixed-point 集合，比较 stable / strict BW / current BU / generalized BU
- 检查 current BU 与导数参考、尾部归一化诊断之间的关系
- 在同一点比较 `mu_s/mu_q=0.2, 0.55`，并检查 charged K/pi 在
  `mu_K = 0` 与 signed flavor difference 规则下的差异

说明：
- 这是分析脚本，不修改正式 workflow 或 kernel
- 输出目录默认写入 `data/outputs/results/relaxtime/analysis/bu_meson_density_literature_audit/`
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "workflow_apps", "WorkflowParamAdapters.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .Models: solve_gap_and_meson_point,
               solve_gap_and_meson_density_point,
               solve_gap_and_strict_bw_meson_density_point,
               solve_gap_and_phase_shift_meson_density_point,
               solve_phase_shift_meson_density_from_meson_point,
               solve_phase_shift_derivative_reference_from_meson_point
using .WorkflowParamAdapters: normalize_quark_params, normalize_thermo_params
using Main.GaussLegendre: gauleg
using Main.MesonDensity: bose_distribution, meson_degeneracy
using Main.MesonMass: ensure_quark_params_has_A
using Main.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using Main.PolarizationAniso: polarization_with_width
using Main.MesonPropagator: meson_propagator_simple

const DEFAULT_OUTDIR = joinpath(
    PROJECT_ROOT,
    "data", "outputs", "results", "relaxtime", "analysis", "bu_meson_density_literature_audit",
)

const NEUTRAL_POINTS = [
    (label="muB0_T208", T_MeV=208.0, muq_MeV=0.0),
    (label="muB0_T220", T_MeV=220.0, muq_MeV=0.0),
]

const CHARGED_POINT = (
    label="charged_diag_T170_muq80",
    T_MeV=170.0,
    muq_MeV=80.0,
    mu_pi_MeV=100.0,
)

const CHARGED_FLAVOR_PROFILES = (
    "bu2020_mu_s_0p2",
    "friesen2019_mu_s_0p55",
)

const STABLE_Q_NODES = 256
const STRICT_BW_QMAX = 12.0
const STRICT_BW_STAGE1_Q_NODES = 48
const STRICT_BW_STAGE2_Q_NODES = 8
const STRICT_BW_OMEGA_MAX = 10.0
const STRICT_BW_STAGE1_OMEGA_NODES = 48
const STRICT_BW_STAGE2_OMEGA_NODES = 16
const PHASE_QMAX = 12.0
const PHASE_Q_NODES = 48
const PHASE_OMEGA_MIN = 0.05
const PHASE_OMEGA_MAX = 10.0
const PHASE_OMEGA_NODES = 48
const PHASE_ETA = 1e-6
const CHARGED_BOSE_X_MIN = 0.05

@inline function _fmt(x)
    x isa Bool && return x ? "true" : "false"
    x isa Symbol && return String(x)
    x isa AbstractString && return replace(replace(x, ',' => ';'), '\n' => ' ')
    return string(x)
end

@inline _env_float(name::AbstractString, default::Real) = parse(Float64, get(ENV, name, string(default)))
@inline _env_int(name::AbstractString, default::Integer) = parse(Int, get(ENV, name, string(default)))

function _charged_phase_settings()
    return (
        qmax=_env_float("BU_AUDIT_QMAX", PHASE_QMAX),
        q_nodes=_env_int("BU_AUDIT_Q_NODES", PHASE_Q_NODES),
        omega_min=_env_float("BU_AUDIT_OMEGA_MIN", PHASE_OMEGA_MIN),
        omega_max=_env_float("BU_AUDIT_OMEGA_MAX", PHASE_OMEGA_MAX),
        omega_nodes=_env_int("BU_AUDIT_OMEGA_NODES", PHASE_OMEGA_NODES),
        eta=_env_float("BU_AUDIT_ETA", PHASE_ETA),
        bose_x_min=_env_float("BU_AUDIT_BOSE_X_MIN", CHARGED_BOSE_X_MIN),
    )
end

@inline function _complex_phase(z::ComplexF64)::Float64
    return atan(imag(z), real(z))
end

function _unwrap_phases(phases::Vector{Float64})
    out = similar(phases)
    isempty(phases) && return out
    out[1] = phases[1]
    shift = 0.0
    for i in 2:length(phases)
        Δ = phases[i] - phases[i - 1]
        if Δ > π
            shift -= 2π
        elseif Δ < -π
            shift += 2π
        end
        out[i] = phases[i] + shift
    end
    return out
end

@inline _gbu_phase_function(δ::Float64) = δ - 0.5 * sin(2.0 * δ)

function _simple_meson_pol_params(meson::Symbol, qp)
    if meson === :pi
        return (
            channel=:P,
            m1=Float64(qp.m.u), m2=Float64(qp.m.u),
            μ1=Float64(qp.μ.u), μ2=Float64(qp.μ.u),
            A1=Float64(qp.A.u), A2=Float64(qp.A.u),
            num_s_quark=0,
        )
    elseif meson === :K
        return (
            channel=:P,
            m1=Float64(qp.m.u), m2=Float64(qp.m.s),
            μ1=Float64(qp.μ.u), μ2=Float64(qp.μ.s),
            A1=Float64(qp.A.u), A2=Float64(qp.A.s),
            num_s_quark=1,
        )
    end
    throw(ArgumentError("Unsupported diagnostic meson: $(meson)"))
end

function _build_k_coeffs(qp)
    G_u = calculate_G_from_A(Float64(qp.A.u), Float64(qp.m.u))
    G_s = calculate_G_from_A(Float64(qp.A.s), Float64(qp.m.s))
    return calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
end

function _propagator_phase(meson::Symbol, ω::Float64, q::Float64, qp, tp, K_coeffs; eta::Float64)
    pol = _simple_meson_pol_params(meson, qp)
    Π_re, Π_im = polarization_with_width(
        pol.channel, ω, 2.0 * eta, q,
        pol.m1, pol.m2,
        pol.μ1, pol.μ2,
        Float64(tp.T), Float64(tp.Φ), Float64(tp.Φbar), Float64(tp.ξ),
        pol.A1, pol.A2, pol.num_s_quark,
    )
    Π = ComplexF64(Π_re, Π_im)
    D = meson_propagator_simple(meson, K_coeffs, Π)
    return _complex_phase(D)
end

function _omega_integral_for_phase(phase_values::Vector{Float64}, omega_grid, omega_w, T::Float64)
    total = 0.0
    @inbounds for iω in eachindex(omega_grid, omega_w, phase_values)
        gω = bose_distribution(Float64(omega_grid[iω]), 0.0, T)
        total += omega_w[iω] * gω * (1.0 + gω) * phase_values[iω]
    end
    return total / (2.0 * π * T)
end

function _tail_normalized_density(meson::Symbol, qp, tp;
                                  qmax::Float64=PHASE_QMAX,
                                  q_nodes::Int=PHASE_Q_NODES,
                                  omega_min::Float64=PHASE_OMEGA_MIN,
                                  omega_max::Float64=PHASE_OMEGA_MAX,
                                  omega_nodes::Int=PHASE_OMEGA_NODES,
                                  eta::Float64=PHASE_ETA)
    K_coeffs = _build_k_coeffs(qp)
    q_grid, q_w = gauleg(0.0, qmax, q_nodes)
    omega_grid, omega_w = gauleg(omega_min, omega_max, omega_nodes)

    raw_sum = 0.0
    tailnorm_sum = 0.0
    gbu_sum = 0.0
    gbu_tailnorm_sum = 0.0

    @inbounds for iq in eachindex(q_grid, q_w)
        q = q_grid[iq]
        phases = [_propagator_phase(meson, ω, q, qp, tp, K_coeffs; eta=eta) for ω in omega_grid]
        delta = _unwrap_phases(phases)
        delta_tail = delta .- last(delta)
        q_pref = q^2 / (2.0 * π^2)
        raw_sum += q_w[iq] * q_pref * _omega_integral_for_phase(delta, omega_grid, omega_w, Float64(tp.T))
        tailnorm_sum += q_w[iq] * q_pref * _omega_integral_for_phase(delta_tail, omega_grid, omega_w, Float64(tp.T))
        gbu_sum += q_w[iq] * q_pref * _omega_integral_for_phase(_gbu_phase_function.(delta), omega_grid, omega_w, Float64(tp.T))
        gbu_tailnorm_sum += q_w[iq] * q_pref * _omega_integral_for_phase(_gbu_phase_function.(delta_tail), omega_grid, omega_w, Float64(tp.T))
    end

    pref = Float64(meson_degeneracy(meson))
    return (
        raw=pref * raw_sum,
        tailnorm=pref * tailnorm_sum,
        gbu_raw=pref * gbu_sum,
        gbu_tailnorm=pref * gbu_tailnorm_sum,
    )
end

function _neutral_point_summary(pt, continuation_state)
    T_fm = pt.T_MeV / ħc_MeV_fm
    muq_fm = pt.muq_MeV / ħc_MeV_fm
    common_kwargs = (
        xi=0.0,
        mesons=(:pi, :K),
        continuation_state=continuation_state,
        mixed_branch_align=:strict_sign_binding,
        p_num=8,
        t_num=4,
        solver_kwargs=(; iterations=20),
        mass_kwargs=(; iterations=20),
    )

    stable = solve_gap_and_meson_density_point(
        T_fm,
        muq_fm;
        common_kwargs...,
        density_kwargs=(; num_q_nodes=STABLE_Q_NODES),
    )
    strict1 = solve_gap_and_strict_bw_meson_density_point(
        T_fm,
        muq_fm;
        common_kwargs...,
        density_kwargs=(;
            stage=:stage1_reduced,
            qmax=STRICT_BW_QMAX,
            q_nodes=STRICT_BW_STAGE1_Q_NODES,
            omega_max=STRICT_BW_OMEGA_MAX,
            omega_nodes=STRICT_BW_STAGE1_OMEGA_NODES,
        ),
    )
    strict2 = solve_gap_and_strict_bw_meson_density_point(
        T_fm,
        muq_fm;
        common_kwargs...,
        density_kwargs=(;
            stage=:stage2_qpole,
            qmax=STRICT_BW_QMAX,
            q_nodes=STRICT_BW_STAGE2_Q_NODES,
            omega_max=STRICT_BW_OMEGA_MAX,
            omega_nodes=STRICT_BW_STAGE2_OMEGA_NODES,
            solver_iterations=12,
            pole_residual_norm_max=1e-4,
            pole_require_converged=true,
        ),
    )
    current = solve_gap_and_phase_shift_meson_density_point(
        T_fm,
        muq_fm;
        common_kwargs...,
        density_kwargs=(;
            scheme=:current,
            qmax=PHASE_QMAX,
            q_nodes=PHASE_Q_NODES,
            omega_min=PHASE_OMEGA_MIN,
            omega_max=PHASE_OMEGA_MAX,
            omega_nodes=PHASE_OMEGA_NODES,
            eta=PHASE_ETA,
        ),
    )
    gbu = solve_gap_and_phase_shift_meson_density_point(
        T_fm,
        muq_fm;
        common_kwargs...,
        density_kwargs=(;
            scheme=:gbu_reference,
            qmax=PHASE_QMAX,
            q_nodes=PHASE_Q_NODES,
            omega_min=PHASE_OMEGA_MIN,
            omega_max=PHASE_OMEGA_MAX,
            omega_nodes=PHASE_OMEGA_NODES,
            eta=PHASE_ETA,
        ),
    )
    derivative_current = solve_phase_shift_derivative_reference_from_meson_point(
        current;
        scheme=:current,
        qmax=PHASE_QMAX,
        q_nodes=PHASE_Q_NODES,
        omega_min=PHASE_OMEGA_MIN,
        omega_max=PHASE_OMEGA_MAX,
        omega_nodes=PHASE_OMEGA_NODES,
        eta=PHASE_ETA,
    )
    derivative_gbu = solve_phase_shift_derivative_reference_from_meson_point(
        current;
        scheme=:gbu_reference,
        qmax=PHASE_QMAX,
        q_nodes=PHASE_Q_NODES,
        omega_min=PHASE_OMEGA_MIN,
        omega_max=PHASE_OMEGA_MAX,
        omega_nodes=PHASE_OMEGA_NODES,
        eta=PHASE_ETA,
    )

    tp = normalize_thermo_params(current.thermo_params)
    qp = ensure_quark_params_has_A(normalize_quark_params(current.quark_params), tp)
    pi_diag = _tail_normalized_density(:pi, qp, tp)
    k_diag = _tail_normalized_density(:K, qp, tp)

    return (
        next_continuation_state=current.continuation_state,
        rows_regime=[
            (point=pt.label, regime="stable", n_pi=stable.meson_density.n_pi, n_K=stable.meson_density.n_K, kpi_ratio=stable.meson_density.kpi_ratio),
            (point=pt.label, regime="strict_bw_stage1", n_pi=strict1.strict_bw_meson_density.n_pi, n_K=strict1.strict_bw_meson_density.n_K, kpi_ratio=strict1.strict_bw_meson_density.kpi_ratio),
            (point=pt.label, regime="strict_bw_stage2", n_pi=strict2.strict_bw_meson_density.n_pi, n_K=strict2.strict_bw_meson_density.n_K, kpi_ratio=strict2.strict_bw_meson_density.kpi_ratio),
            (point=pt.label, regime="phase_shift_current", n_pi=current.phase_shift_meson_density.n_pi, n_K=current.phase_shift_meson_density.n_K, kpi_ratio=current.phase_shift_meson_density.kpi_ratio),
            (point=pt.label, regime="phase_shift_gbu", n_pi=gbu.phase_shift_meson_density.n_pi, n_K=gbu.phase_shift_meson_density.n_K, kpi_ratio=gbu.phase_shift_meson_density.kpi_ratio),
            (point=pt.label, regime="phase_shift_current_derivative_ref", n_pi=derivative_current.n_pi, n_K=derivative_current.n_K, kpi_ratio=derivative_current.kpi_ratio),
            (point=pt.label, regime="phase_shift_gbu_derivative_ref", n_pi=derivative_gbu.n_pi, n_K=derivative_gbu.n_K, kpi_ratio=derivative_gbu.kpi_ratio),
        ],
        rows_diag=[
            (point=pt.label, meson="pi", density_raw=pi_diag.raw, density_tailnorm=pi_diag.tailnorm, density_gbu_raw=pi_diag.gbu_raw, density_gbu_tailnorm=pi_diag.gbu_tailnorm),
            (point=pt.label, meson="K", density_raw=k_diag.raw, density_tailnorm=k_diag.tailnorm, density_gbu_raw=k_diag.gbu_raw, density_gbu_tailnorm=k_diag.gbu_tailnorm),
        ],
    )
end

function _charged_muK_rows()
    T_fm = CHARGED_POINT.T_MeV / ħc_MeV_fm
    muq_fm = CHARGED_POINT.muq_MeV / ħc_MeV_fm
    mu_pi_fm = CHARGED_POINT.mu_pi_MeV / ħc_MeV_fm
    settings = _charged_phase_settings()
    rows = NamedTuple[]

    for profile_name in CHARGED_FLAVOR_PROFILES
        profile = Models.FlavorChemicalProfiles.load_flavor_chemical_profile(profile=profile_name)
        flavor_mev = Models.FlavorChemicalProfiles.flavor_mu_profile_MeV(profile, CHARGED_POINT.muq_MeV)
        flavor_fm = Models.FlavorChemicalProfiles.flavor_mu_profile_fm(profile, CHARGED_POINT.muq_MeV)
        muK_plus_signed_fm = flavor_fm.mu_u_fm - flavor_fm.mu_s_fm
        muK_minus_signed_fm = -muK_plus_signed_fm
        equilibrium_group = "$(CHARGED_POINT.label)__$(profile.profile_name)"

        # One and only one upstream state solve per flavor profile. Both charge
        # channels and both mu_K rules below are post-processing of this state.
        meson_point = solve_gap_and_meson_point(
            T_fm,
            muq_fm;
            xi=0.0,
            mesons=(:pi_plus, :K_plus, :pi_minus, :K_minus),
            mixed_branch_align=:strict_sign_binding,
            flavor_mu_override=(flavor_fm.mu_u_fm, flavor_fm.mu_d_fm, flavor_fm.mu_s_fm),
            p_num=8,
            t_num=4,
            solver_kwargs=(; iterations=20),
            mass_kwargs=(; iterations=20),
        )

        function _charged_summary(pi_channel::Symbol, k_channel::Symbol, mu_K_fm::Float64)
            return solve_phase_shift_meson_density_from_meson_point(
                meson_point;
                pi_channel=pi_channel,
                k_channel=k_channel,
                μ_pi=mu_pi_fm,
                μ_K=mu_K_fm,
                d_pi=1,
                d_K=1,
                scheme=:gbu_reference,
                qmax=settings.qmax,
                q_nodes=settings.q_nodes,
                omega_min=settings.omega_min,
                omega_max=settings.omega_max,
                omega_nodes=settings.omega_nodes,
                eta=settings.eta,
                density_policy=:x_min_cut,
                bose_x_min=settings.bose_x_min,
            )
        end

        plus_mu0 = _charged_summary(:pi_plus, :K_plus, 0.0)
        minus_mu0 = _charged_summary(:pi_minus, :K_minus, 0.0)
        plus_signed = _charged_summary(:pi_plus, :K_plus, muK_plus_signed_fm)
        minus_signed = _charged_summary(:pi_minus, :K_minus, muK_minus_signed_fm)

        function _row(channel::String, mu_rule::String, mu_formula::String,
                      pi_channel::Symbol, k_channel::Symbol,
                      mu_K_fm::Float64, density, peer_density)
            pi_mass_fm = Float64(meson_point.meson_results[pi_channel].mass)
            k_mass_fm = Float64(meson_point.meson_results[k_channel].mass)
            is_plus = channel == "Kplus_over_piplus"
            n_pi_plus = is_plus ? density.n_pi : peer_density.n_pi
            n_K_plus = is_plus ? density.n_K : peer_density.n_K
            n_pi_minus = is_plus ? peer_density.n_pi : density.n_pi
            n_K_minus = is_plus ? peer_density.n_K : density.n_K
            strict_status = settings.omega_min > max(mu_pi_fm, mu_K_fm) ? :strict_safe : :strict_unsafe_support
            return (
                point=CHARGED_POINT.label,
                equilibrium_group=equilibrium_group,
                upstream_state_id=equilibrium_group,
                scenario_role="literature_sensitivity_not_conserved_charge_solution",
                flavor_profile=profile.profile_name,
                flavor_source_tag=profile.source_tag,
                flavor_apply_to_equilibrium=profile.apply_to_equilibrium,
                mu_s_over_muq=profile.mu_s_over_muq,
                T_MeV=CHARGED_POINT.T_MeV,
                muq_MeV=CHARGED_POINT.muq_MeV,
                muB_MeV=3.0 * CHARGED_POINT.muq_MeV,
                mu_u_MeV=flavor_mev.mu_u_MeV,
                mu_d_MeV=flavor_mev.mu_d_MeV,
                mu_s_MeV=flavor_mev.mu_s_MeV,
                channel=channel,
                pi_channel=String(pi_channel),
                k_channel=String(k_channel),
                mu_rule=mu_rule,
                mu_formula=mu_formula,
                mu_pi_MeV=CHARGED_POINT.mu_pi_MeV,
                mu_K_MeV=mu_K_fm * ħc_MeV_fm,
                m_pi_MeV=pi_mass_fm * ħc_MeV_fm,
                m_K_MeV=k_mass_fm * ħc_MeV_fm,
                m_pi_minus_mu_pi_MeV=(pi_mass_fm - mu_pi_fm) * ħc_MeV_fm,
                m_K_minus_mu_K_MeV=(k_mass_fm - mu_K_fm) * ħc_MeV_fm,
                n_pi=density.n_pi,
                n_K=density.n_K,
                kpi_ratio=density.kpi_ratio,
                n_pi_plus=n_pi_plus,
                n_K_plus=n_K_plus,
                n_pi_minus=n_pi_minus,
                n_K_minus=n_K_minus,
                channel_status=density.status,
                pi_status=density.pi_density.status,
                k_status=density.k_density.status,
                strict_requested_window_status=strict_status,
                unsafe_bose_count=density.unsafe_bose_count,
                min_E_minus_mu=density.min_E_minus_mu,
                density_policy=density.density_policy,
                bose_x_min=density.bose_x_min,
                omega_min_effective_pi=density.pi_density.omega_min_effective,
                omega_min_effective_K=density.k_density.omega_min_effective,
                message=density.message,
                qmax=settings.qmax,
                q_nodes=settings.q_nodes,
                omega_min=settings.omega_min,
                omega_max=settings.omega_max,
                omega_nodes=settings.omega_nodes,
                eta=settings.eta,
            )
        end

        append!(rows, (
            _row("Kplus_over_piplus", "muK_zero", "zero", :pi_plus, :K_plus, 0.0, plus_mu0, minus_mu0),
            _row("Kminus_over_piminus", "muK_zero", "zero", :pi_minus, :K_minus, 0.0, minus_mu0, plus_mu0),
            _row("Kplus_over_piplus", "signed_flavor_difference", "mu_u_minus_mu_s", :pi_plus, :K_plus, muK_plus_signed_fm, plus_signed, minus_signed),
            _row("Kminus_over_piminus", "signed_flavor_difference", "mu_s_minus_mu_u", :pi_minus, :K_minus, muK_minus_signed_fm, minus_signed, plus_signed),
        ))
    end

    return rows
end

function _write_csv(path::String, header::Vector{String}, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(header, ','))
        for row in rows
            println(io, join((_fmt(getproperty(row, Symbol(col))) for col in header), ','))
        end
    end
end

function main()
    outdir = isempty(ARGS) ? DEFAULT_OUTDIR : abspath(ARGS[1])
    mkpath(outdir)

    regime_rows = NamedTuple[]
    diag_rows = NamedTuple[]
    continuation_state = nothing
    run_neutral = lowercase(get(ENV, "BU_AUDIT_RUN_NEUTRAL", "true")) in ("1", "true", "yes")
    if run_neutral
        for pt in NEUTRAL_POINTS
            result = _neutral_point_summary(pt, continuation_state)
            append!(regime_rows, result.rows_regime)
            append!(diag_rows, result.rows_diag)
            continuation_state = result.next_continuation_state
        end
    end

    charged_rows = _charged_muK_rows()

    _write_csv(
        joinpath(outdir, "fixedpoint_regime_summary.csv"),
        ["point", "regime", "n_pi", "n_K", "kpi_ratio"],
        regime_rows,
    )
    _write_csv(
        joinpath(outdir, "fixedpoint_phase_normalization_diag.csv"),
        ["point", "meson", "density_raw", "density_tailnorm", "density_gbu_raw", "density_gbu_tailnorm"],
        diag_rows,
    )
    _write_csv(
        joinpath(outdir, "charged_muK_rule_sensitivity.csv"),
        [
            "point", "equilibrium_group", "upstream_state_id", "scenario_role",
            "flavor_profile", "flavor_source_tag", "flavor_apply_to_equilibrium", "mu_s_over_muq",
            "T_MeV", "muq_MeV", "muB_MeV", "mu_u_MeV", "mu_d_MeV", "mu_s_MeV",
            "channel", "pi_channel", "k_channel", "mu_rule", "mu_formula", "mu_pi_MeV", "mu_K_MeV",
            "m_pi_MeV", "m_K_MeV", "m_pi_minus_mu_pi_MeV", "m_K_minus_mu_K_MeV",
            "n_pi", "n_K", "kpi_ratio", "n_pi_plus", "n_K_plus", "n_pi_minus", "n_K_minus",
            "channel_status", "pi_status", "k_status", "strict_requested_window_status",
            "unsafe_bose_count", "min_E_minus_mu", "density_policy", "bose_x_min",
            "omega_min_effective_pi", "omega_min_effective_K", "message",
            "qmax", "q_nodes", "omega_min", "omega_max", "omega_nodes", "eta",
        ],
        charged_rows,
    )

    println("Wrote audit outputs to: ", outdir)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
