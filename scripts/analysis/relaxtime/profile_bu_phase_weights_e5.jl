"""
Phase E5: 剖面化输出 BU 相位与相位函数权重。

目标：
- 在少量 `(T, q)` 取样点上，直接输出：
  - `delta(omega, q)`
  - `F_current(delta) = delta`
  - `F_gbu(delta) = delta - 0.5 sin(2 delta)`
  - 对应的被积函数权重
- 配合 threshold below / above 标记，观察 generalized BU 究竟压掉了哪类结构

说明：
- 该脚本是分析层工具，不改正式 helper / workflow
- 输出 detail CSV + summary CSV
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
    "bu_phase_profile_e5",
)

const DEFAULT_TS_MEV = [208.0, 210.0, 212.0]
const DEFAULT_Q_VALUES = [0.0, 2.0, 4.0]
const DEFAULT_OMEGA_MIN = 0.05
const DEFAULT_OMEGA_MAX = 10.0
const DEFAULT_OMEGA_STEP = 0.05
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

    detail_path = joinpath(outdir, "bu_phase_profile_detail.csv")
    summary_path = joinpath(outdir, "bu_phase_profile_summary.csv")

    continuation_state = nothing
    summary_rows = NamedTuple[]

    open(detail_path, "w") do io
        println(io, "T_MeV,meson,q,omega,omega_thr,region,mass_ref,gamma_ref,threshold_ref,delta,F_current,F_gbu,bose_kernel,integrand_current,integrand_gbu")

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
                    raw_phases = [_propagator_phase(meson, ω, q, qp, tp, K_coeffs; eta=eta) for ω in omegas]
                    delta = _unwrap_phases(raw_phases)

                    current_below = 0.0
                    current_above = 0.0
                    gbu_below = 0.0
                    gbu_above = 0.0
                    delta_max = maximum(delta)
                    delta_min = minimum(delta)
                    omega_thr = hypot(q, Float64(mr.threshold))

                    for (ω, δ) in zip(omegas, delta)
                        F_current = δ
                        F_gbu = _gbu_phase_function(δ)
                        kernel = bose_distribution(Float64(ω), 0.0, Float64(tp.T)) *
                                 (1.0 + bose_distribution(Float64(ω), 0.0, Float64(tp.T))) /
                                 (2.0 * π * Float64(tp.T))
                        ic = kernel * F_current
                        ig = kernel * F_gbu
                        is_above = ω > omega_thr
                        region = is_above ? "above_threshold" : "below_threshold"
                        if is_above
                            current_above += ic
                            gbu_above += ig
                        else
                            current_below += ic
                            gbu_below += ig
                        end
                        println(io, join((
                            _fmt(T_MeV), String(meson), _fmt(q), _fmt(ω), _fmt(omega_thr), region,
                            _fmt(Float64(mr.mass)), _fmt(Float64(mr.gamma)), _fmt(Float64(mr.threshold)),
                            _fmt(δ), _fmt(F_current), _fmt(F_gbu), _fmt(kernel), _fmt(ic), _fmt(ig)
                        ), ','))
                    end

                    push!(summary_rows, (
                        T_MeV=T_MeV,
                        meson=String(meson),
                        q=q,
                        omega_thr=omega_thr,
                        mass_ref=Float64(mr.mass),
                        gamma_ref=Float64(mr.gamma),
                        threshold_ref=Float64(mr.threshold),
                        delta_min=delta_min,
                        delta_max=delta_max,
                        current_below=current_below,
                        current_above=current_above,
                        gbu_below=gbu_below,
                        gbu_above=gbu_above,
                    ))
                end
            end
        end
    end

    open(summary_path, "w") do io
        println(io, "T_MeV,meson,q,omega_thr,mass_ref,gamma_ref,threshold_ref,delta_min,delta_max,current_below,current_above,gbu_below,gbu_above")
        for r in summary_rows
            println(io, join((_fmt(getproperty(r, k)) for k in (
                :T_MeV, :meson, :q, :omega_thr, :mass_ref, :gamma_ref, :threshold_ref,
                :delta_min, :delta_max, :current_below, :current_above, :gbu_below, :gbu_above
            )), ','))
        end
    end

    println("Wrote Phase E5 profile detail CSV: ", detail_path)
    println("Wrote Phase E5 profile summary CSV: ", summary_path)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
