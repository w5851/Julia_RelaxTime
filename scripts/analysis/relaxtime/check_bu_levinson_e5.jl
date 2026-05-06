"""
Phase E5: Levinson / phase normalization 诊断。

目标：
- 固定少量 `(T, q)` 采样点
- 检查 `delta(omega, q)` 的：
  1. 低能端基线
  2. 高能端尾部
  3. 总相位变化量
  4. 与整数倍 `pi` 的接近程度
- 为判断当前 `arg D` 总相位是否已足够接近“更严格 scattering phase shift”提供门禁证据
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

const DEFAULT_OUTPUT = joinpath(
    PROJECT_ROOT,
    "data",
    "outputs",
    "results",
    "relaxtime",
    "analysis",
    "meson_density_phase_e5_levinson_check.csv",
)

const DEFAULT_TS_MEV = [208.0, 210.0, 212.0]
const DEFAULT_Q_VALUES = [0.0, 2.0, 4.0]
const DEFAULT_OMEGA_MIN = 0.05
const DEFAULT_OMEGA_MAX = 10.0
const DEFAULT_OMEGA_STEP = 0.02
const DEFAULT_ETA = 1e-6

@inline function _fmt(x)
    x isa Bool && return (x ? "true" : "false")
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

function _selected_temperatures()
    raw = strip(get(ENV, "MESON_E5_T_VALUES", ""))
    isempty(raw) && return DEFAULT_TS_MEV
    return _parse_float_list(raw)
end

function _selected_q_values()
    raw = strip(get(ENV, "MESON_E5_Q_VALUES", ""))
    isempty(raw) && return DEFAULT_Q_VALUES
    return _parse_float_list(raw)
end

function _omega_grid()
    ωmin = parse(Float64, get(ENV, "MESON_E5_OMEGA_MIN", string(DEFAULT_OMEGA_MIN)))
    ωmax = parse(Float64, get(ENV, "MESON_E5_OMEGA_MAX", string(DEFAULT_OMEGA_MAX)))
    ωstep = parse(Float64, get(ENV, "MESON_E5_OMEGA_STEP", string(DEFAULT_OMEGA_STEP)))
    ωstep > 0 || throw(ArgumentError("MESON_E5_OMEGA_STEP must be positive"))
    ωmax > ωmin || throw(ArgumentError("MESON_E5_OMEGA_MAX must exceed MESON_E5_OMEGA_MIN"))
    collect(ωmin:ωstep:ωmax)
end

@inline function _complex_phase(z::ComplexF64)::Float64
    atan(imag(z), real(z))
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

function _build_k_coeffs(qp)
    G_u = calculate_G_from_A(Float64(qp.A.u), Float64(qp.m.u))
    G_s = calculate_G_from_A(Float64(qp.A.s), Float64(qp.m.s))
    calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
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
    _complex_phase(D)
end

@inline function _nearest_pi_multiple(x::Float64)
    n = round(Int, x / π)
    return n, n * π, abs(x - n * π)
end

function main()
    output = isempty(ARGS) ? DEFAULT_OUTPUT : abspath(ARGS[1])
    mkpath(dirname(output))

    Ts_MeV = _selected_temperatures()
    qs = _selected_q_values()
    omegas = _omega_grid()
    eta = parse(Float64, get(ENV, "MESON_E5_ETA", string(DEFAULT_ETA)))

    continuation_state = nothing
    rows = NamedTuple[]

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

        tp = normalize_thermo_params(res.thermo_params)
        qp = ensure_quark_params_has_A(normalize_quark_params(res.quark_params), tp)
        K_coeffs = _build_k_coeffs(qp)

        for meson in (:pi, :K)
            mr = res.meson_results[meson]
            for q in qs
                raw = [_propagator_phase(meson, ω, q, qp, tp, K_coeffs; eta=eta) for ω in omegas]
                delta = _unwrap_phases(raw)
                δ_lo = first(delta)
                δ_hi = last(delta)
                δ_span = δ_hi - δ_lo
                n_lo, npi_lo, dist_lo = _nearest_pi_multiple(δ_lo)
                n_hi, npi_hi, dist_hi = _nearest_pi_multiple(δ_hi)
                n_span, npi_span, dist_span = _nearest_pi_multiple(δ_span)
                omega_thr = hypot(q, Float64(mr.threshold))

                push!(rows, (
                    T_MeV=T_MeV,
                    meson=String(meson),
                    q=q,
                    mass_ref=Float64(mr.mass),
                    gamma_ref=Float64(mr.gamma),
                    threshold_ref=Float64(mr.threshold),
                    omega_thr=omega_thr,
                    delta_low=δ_lo,
                    delta_high=δ_hi,
                    delta_span=δ_span,
                    low_pi_multiple=n_lo,
                    low_pi_target=npi_lo,
                    low_pi_distance=dist_lo,
                    high_pi_multiple=n_hi,
                    high_pi_target=npi_hi,
                    high_pi_distance=dist_hi,
                    span_pi_multiple=n_span,
                    span_pi_target=npi_span,
                    span_pi_distance=dist_span,
                    delta_min=minimum(delta),
                    delta_max=maximum(delta),
                ))
            end
        end
    end

    open(output, "w") do io
        println(io, "T_MeV,meson,q,mass_ref,gamma_ref,threshold_ref,omega_thr,delta_low,delta_high,delta_span,low_pi_multiple,low_pi_target,low_pi_distance,high_pi_multiple,high_pi_target,high_pi_distance,span_pi_multiple,span_pi_target,span_pi_distance,delta_min,delta_max")
        for r in rows
            println(io, join((_fmt(getproperty(r, k)) for k in (
                :T_MeV, :meson, :q, :mass_ref, :gamma_ref, :threshold_ref, :omega_thr,
                :delta_low, :delta_high, :delta_span,
                :low_pi_multiple, :low_pi_target, :low_pi_distance,
                :high_pi_multiple, :high_pi_target, :high_pi_distance,
                :span_pi_multiple, :span_pi_target, :span_pi_distance,
                :delta_min, :delta_max
            )), ','))
        end
    end

    println("Wrote Phase E5 Levinson CSV: ", output)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
