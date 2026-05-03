"""
Phase E5: 围绕极点与阈值的局部窗口剖面。

目标：
- 固定少量 `(T, q)` 采样点
- 分别围绕
  1. pole window: `omega ~ mass_ref`
  2. threshold window: `omega ~ omega_thr(q)`
  输出局部相位与被积函数
- 判断 generalized BU 主要压掉的是 pole-neighborhood 还是 threshold-neighborhood 的结构
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "workflow_apps", "WorkflowParamAdapters.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .Models: solve_gap_and_meson_point
using .WorkflowParamAdapters: normalize_quark_params, normalize_thermo_params
using Main.MesonDensity: bose_distribution
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
    "bu_phase_local_windows_e5",
)

const DEFAULT_TS_MEV = [208.0, 210.0, 212.0]
const DEFAULT_Q_VALUES = [0.0, 2.0]
const DEFAULT_OMEGA_MIN = 0.05
const DEFAULT_OMEGA_MAX = 10.0
const DEFAULT_OMEGA_STEP = 0.02
const DEFAULT_WINDOW_HALF_WIDTH = 0.4
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

@inline _gbu_phase_function(δ::Float64) = δ - 0.5 * sin(2.0 * δ)

function main()
    outdir = isempty(ARGS) ? DEFAULT_OUTDIR : abspath(ARGS[1])
    mkpath(outdir)

    Ts_MeV = _selected_temperatures()
    qs = _selected_q_values()
    omegas = _omega_grid()
    eta = parse(Float64, get(ENV, "MESON_E5_ETA", string(DEFAULT_ETA)))
    half_width = parse(Float64, get(ENV, "MESON_E5_WINDOW_HALF_WIDTH", string(DEFAULT_WINDOW_HALF_WIDTH)))

    detail_path = joinpath(outdir, "bu_phase_local_windows_detail.csv")
    summary_path = joinpath(outdir, "bu_phase_local_windows_summary.csv")

    continuation_state = nothing
    summary_rows = NamedTuple[]

    open(detail_path, "w") do io
        println(io, "T_MeV,meson,q,window_type,omega,window_center,omega_thr,mass_ref,gamma_ref,threshold_ref,delta,F_current,F_gbu,integrand_current,integrand_gbu")

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
                    omega_thr = hypot(q, Float64(mr.threshold))
                    pole_center = Float64(mr.mass)
                    for (window_type, center) in (("pole_window", pole_center), ("threshold_window", omega_thr))
                        lo = center - half_width
                        hi = center + half_width
                        idx = findall(ω -> (ω >= lo && ω <= hi), omegas)
                        isempty(idx) && continue
                        ωsel = omegas[idx]
                        raw = [_propagator_phase(meson, ω, q, qp, tp, K_coeffs; eta=eta) for ω in ωsel]
                        delta = _unwrap_phases(raw)

                        current_sum = 0.0
                        gbu_sum = 0.0
                        for (ω, δ) in zip(ωsel, delta)
                            F_current = δ
                            F_gbu = _gbu_phase_function(δ)
                            gω = bose_distribution(Float64(ω), 0.0, Float64(tp.T))
                            ic = gω * (1.0 + gω) * F_current / (2.0 * π * Float64(tp.T))
                            ig = gω * (1.0 + gω) * F_gbu / (2.0 * π * Float64(tp.T))
                            current_sum += ic
                            gbu_sum += ig
                            println(io, join((
                                _fmt(T_MeV), String(meson), _fmt(q), window_type, _fmt(ω), _fmt(center), _fmt(omega_thr),
                                _fmt(Float64(mr.mass)), _fmt(Float64(mr.gamma)), _fmt(Float64(mr.threshold)),
                                _fmt(δ), _fmt(F_current), _fmt(F_gbu), _fmt(ic), _fmt(ig)
                            ), ','))
                        end

                        push!(summary_rows, (
                            T_MeV=T_MeV,
                            meson=String(meson),
                            q=q,
                            window_type=window_type,
                            window_center=center,
                            omega_thr=omega_thr,
                            mass_ref=Float64(mr.mass),
                            gamma_ref=Float64(mr.gamma),
                            threshold_ref=Float64(mr.threshold),
                            delta_min=minimum(delta),
                            delta_max=maximum(delta),
                            current_window_sum=current_sum,
                            gbu_window_sum=gbu_sum,
                        ))
                    end
                end
            end
        end
    end

    open(summary_path, "w") do io
        println(io, "T_MeV,meson,q,window_type,window_center,omega_thr,mass_ref,gamma_ref,threshold_ref,delta_min,delta_max,current_window_sum,gbu_window_sum")
        for r in summary_rows
            println(io, join((_fmt(getproperty(r, k)) for k in (
                :T_MeV, :meson, :q, :window_type, :window_center, :omega_thr,
                :mass_ref, :gamma_ref, :threshold_ref, :delta_min, :delta_max,
                :current_window_sum, :gbu_window_sum
            )), ','))
        end
    end

    println("Wrote Phase E5 local-window detail CSV: ", detail_path)
    println("Wrote Phase E5 local-window summary CSV: ", summary_path)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
