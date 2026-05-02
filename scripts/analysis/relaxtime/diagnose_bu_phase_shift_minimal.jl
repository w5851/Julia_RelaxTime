"""
Phase E1: 最小 BU 相移诊断脚本。

目标：
- 复用现有 meson workflow 获取 `quark_params` / `thermo_params`
- 通过 `PolarizationAniso + MesonPropagator` 链构造简单介子道传播子
- 在少量温点、固定 q 下输出 `D_M(omega)` 与 `delta_M(omega)=arg D_M`
- 为后续正式 BU 数密度积分提供相位结构门禁数据

当前固定范围：
- muB = 0
- xi = 0
- mesons = pi, K
- q = 0
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "workflow_apps", "WorkflowParamAdapters.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .Models: solve_gap_and_meson_point
using .WorkflowParamAdapters: normalize_quark_params, normalize_thermo_params
using Main.MesonMass: ensure_quark_params_has_A
using Main.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using Main.PolarizationAniso: polarization_with_width
using Main.MesonPropagator: meson_propagator_simple

const DEFAULT_OUTDIR = joinpath(
    PROJECT_ROOT,
    "data",
    "outputs",
    "results",
    "relaxtime",
    "analysis",
    "bu_phase_shift_minimal",
)

const DEFAULT_T_VALUES_MEV = [208.0, 210.0, 212.0]
const DEFAULT_OMEGA_MIN = 0.05
const DEFAULT_OMEGA_MAX = 4.0
const DEFAULT_OMEGA_STEP = 0.01
const DEFAULT_Q_VALUES = [0.0]
const DEFAULT_ETA = 1e-6

@inline function _fmt(x)
    x isa Bool && return x ? "true" : "false"
    return string(x)
end

function _parse_float_list(raw::AbstractString)
    vals = Float64[]
    for seg in split(raw, ',')
        s = strip(seg)
        isempty(s) && continue
        push!(vals, parse(Float64, s))
    end
    isempty(vals) && throw(ArgumentError("list cannot be empty"))
    return vals
end

@inline function _exp_omega_over_T(omega::Float64, T::Float64)::Float64
    T > 0.0 || return Inf
    x = omega / T
    x > 700.0 && return Inf
    return exp(x)
end

@inline function _bose_weight(omega::Float64, T::Float64)::Float64
    T > 0.0 || return 0.0
    omega > 0.0 || return Inf
    x = omega / T
    x > 700.0 && return 0.0
    return 1.0 / expm1(x)
end

function _selected_temperatures()
    raw = strip(get(ENV, "PHASE_E_T_VALUES", ""))
    isempty(raw) && return DEFAULT_T_VALUES_MEV
    return _parse_float_list(raw)
end

function _selected_q_values()
    raw = strip(get(ENV, "PHASE_E_Q_VALUES", ""))
    isempty(raw) && return DEFAULT_Q_VALUES
    return _parse_float_list(raw)
end

function _omega_grid()
    ωmin = parse(Float64, get(ENV, "PHASE_E_OMEGA_MIN", string(DEFAULT_OMEGA_MIN)))
    ωmax = parse(Float64, get(ENV, "PHASE_E_OMEGA_MAX", string(DEFAULT_OMEGA_MAX)))
    ωstep = parse(Float64, get(ENV, "PHASE_E_OMEGA_STEP", string(DEFAULT_OMEGA_STEP)))
    ωstep > 0 || throw(ArgumentError("PHASE_E_OMEGA_STEP must be positive"))
    ωmax > ωmin || throw(ArgumentError("PHASE_E_OMEGA_MAX must exceed PHASE_E_OMEGA_MIN"))
    return collect(ωmin:ωstep:ωmax)
end

@inline function _propagator_channel(meson::Symbol)::Symbol
    if meson === :pi || meson === :K
        return :P
    elseif meson === :sigma_pi || meson === :sigma_K
        return :S
    end
    throw(ArgumentError("Unsupported simple meson: $(meson)"))
end

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
    throw(ArgumentError("Unsupported simple meson: $(meson)"))
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

function _build_k_coeffs(qp)
    G_u = calculate_G_from_A(Float64(qp.A.u), Float64(qp.m.u))
    G_s = calculate_G_from_A(Float64(qp.A.s), Float64(qp.m.s))
    return calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
end

function _compute_simple_propagator_row(meson::Symbol, ω::Float64, q::Float64, qp, tp, K_coeffs; eta::Float64)
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
    Kc = if meson === :pi
        K_coeffs.K123_plus
    elseif meson === :K
        K_coeffs.K4567_plus
    else
        throw(ArgumentError("Unsupported simple meson: $(meson)"))
    end
    denom = 1.0 - 4.0 * Kc * Π
    return (
        omega=ω,
        q=q,
        exp_omega_over_T=_exp_omega_over_T(ω, Float64(tp.T)),
        bose_g=_bose_weight(ω, Float64(tp.T)),
        Pi_re=real(Π),
        Pi_im=imag(Π),
        denom_re=real(denom),
        denom_im=imag(denom),
        denom_phase=_complex_phase(denom),
        D_re=real(D),
        D_im=imag(D),
        D_abs=abs(D),
        phase_raw=_complex_phase(D),
    )
end

function main()
    outdir = isempty(ARGS) ? DEFAULT_OUTDIR : abspath(ARGS[1])
    mkpath(outdir)

    Ts_MeV = _selected_temperatures()
    qs = _selected_q_values()
    omegas = _omega_grid()
    eta = parse(Float64, get(ENV, "PHASE_E_ETA", string(DEFAULT_ETA)))

    detail_path = joinpath(outdir, "bu_phase_shift_detail.csv")
    summary_path = joinpath(outdir, "bu_phase_shift_summary.csv")

    continuation_state = nothing
    summary_rows = NamedTuple[]

    open(detail_path, "w") do io
        println(io, "T_MeV,xi,meson,q,omega,mass_ref,gamma_ref,threshold_ref,gap_ref,exp_omega_over_T,bose_g,Pi_re,Pi_im,denom_re,denom_im,denom_phase,D_re,D_im,D_abs,phase_raw,phase_unwrapped,phase_shifted_tail")

        for T_MeV in Ts_MeV
            T_fm = T_MeV / ħc_MeV_fm
            res = solve_gap_and_meson_point(
                T_fm,
                0.0;
                xi=0.0,
                mesons=(:pi, :K),
                continuation_state=continuation_state,
                mixed_branch_align=:strict_sign_binding,
                p_num=8,
                t_num=4,
                solver_kwargs=(; iterations=20),
                mass_kwargs=(; iterations=20),
            )
            continuation_state = res.continuation_state

            qp = ensure_quark_params_has_A(
                normalize_quark_params(res.quark_params),
                normalize_thermo_params(res.thermo_params),
            )
            tp = normalize_thermo_params(res.thermo_params)
            K_coeffs = _build_k_coeffs(qp)

            for meson in (:pi, :K)
                mr = res.meson_results[meson]
                threshold_ref = Float64(mr.threshold)
                gap_ref = Float64(mr.gap)

                for q in qs
                    rows = [_compute_simple_propagator_row(meson, ω, q, qp, tp, K_coeffs; eta=eta) for ω in omegas]
                    phase_unwrapped = _unwrap_phases([r.phase_raw for r in rows])
                    tail_ref = phase_unwrapped[end]

                    for (r, δu) in zip(rows, phase_unwrapped)
                        println(io, join((
                            _fmt(T_MeV), _fmt(Float64(tp.ξ)), String(meson), _fmt(r.q), _fmt(r.omega),
                            _fmt(Float64(mr.mass)), _fmt(Float64(mr.gamma)), _fmt(threshold_ref), _fmt(gap_ref),
                            _fmt(r.exp_omega_over_T), _fmt(r.bose_g),
                            _fmt(r.Pi_re), _fmt(r.Pi_im),
                            _fmt(r.denom_re), _fmt(r.denom_im), _fmt(r.denom_phase),
                            _fmt(r.D_re), _fmt(r.D_im), _fmt(r.D_abs),
                            _fmt(r.phase_raw), _fmt(δu), _fmt(δu - tail_ref),
                        ), ','))
                    end

                    push!(summary_rows, (
                        T_MeV=T_MeV,
                        xi=Float64(tp.ξ),
                        meson=String(meson),
                        q=q,
                        mass_ref=Float64(mr.mass),
                        gamma_ref=Float64(mr.gamma),
                        threshold_ref=threshold_ref,
                        gap_ref=gap_ref,
                        phase_min=minimum(phase_unwrapped),
                        phase_max=maximum(phase_unwrapped),
                        phase_span=maximum(phase_unwrapped) - minimum(phase_unwrapped),
                        tail_phase=tail_ref,
                        shifted_phase_min=minimum(δ - tail_ref for δ in phase_unwrapped),
                        shifted_phase_max=maximum(δ - tail_ref for δ in phase_unwrapped),
                        D_abs_max=maximum(r.D_abs for r in rows),
                    ))
                end
            end
        end
    end

    open(summary_path, "w") do io
        println(io, "T_MeV,xi,meson,q,mass_ref,gamma_ref,threshold_ref,gap_ref,phase_min,phase_max,phase_span,tail_phase,shifted_phase_min,shifted_phase_max,D_abs_max")
        for r in summary_rows
            println(io, join((_fmt(getproperty(r, k)) for k in (
                :T_MeV, :xi, :meson, :q, :mass_ref, :gamma_ref, :threshold_ref, :gap_ref,
                :phase_min, :phase_max, :phase_span, :tail_phase, :shifted_phase_min, :shifted_phase_max, :D_abs_max
            )), ','))
        end
    end

    println("Wrote Phase E detail CSV: ", detail_path)
    println("Wrote Phase E summary CSV: ", summary_path)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
