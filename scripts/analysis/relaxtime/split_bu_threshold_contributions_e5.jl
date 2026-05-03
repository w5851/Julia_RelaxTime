"""
Phase E5: 按阈值拆分 BU 相位积分贡献。

目标：
- 固定与 Phase-E3 / E5 一致的 meson workflow / continuation 口径
- 在同一组 `(T, q, omega)` 网格口径下，将 BU 数密度积分拆成：
  1. threshold below
  2. threshold above
- 同时比较：
  - current: `F(delta) = delta`
  - generalized BU: `F(delta) = delta - 0.5 * sin(2 delta)`

说明：
- 本脚本仍是分析层，不修改正式 workflow / helper
- 阈值采用当前 meson workflow 返回的 `threshold_ref = threshold(q=0)`，
  并在每个 `q` 壳上提升为
    `omega_thr(q) = sqrt(q^2 + threshold_ref^2)`
  作为第一版 continuum 分界
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "workflow_apps", "WorkflowParamAdapters.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .Models: solve_gap_and_meson_point
using .WorkflowParamAdapters: normalize_quark_params, normalize_thermo_params
using Main.GaussLegendre: gauleg
using Main.MesonDensity: bose_distribution, meson_degeneracy
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
    "meson_density_phase_e5_threshold_split.csv",
)

const DEFAULT_TS_MEV = [208.0, 210.0, 212.0]
const DEFAULT_QMAX = 12.0
const DEFAULT_Q_NODES = 48
const DEFAULT_OMEGA_MIN = 0.05
const DEFAULT_OMEGA_MAX = 10.0
const DEFAULT_OMEGA_NODES = 48
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
    isempty(vals) && throw(ArgumentError("temperature list cannot be empty"))
    return vals
end

function _selected_temperatures()
    raw = strip(get(ENV, "MESON_E5_T_VALUES", ""))
    isempty(raw) && return DEFAULT_TS_MEV
    return _parse_float_list(raw)
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

function _omega_phase_integral(meson::Symbol, q::Float64, ωa::Float64, ωb::Float64,
                               qp, tp, K_coeffs; omega_nodes::Int, eta::Float64)
    if !(ωb > ωa)
        return (current=0.0, gbu=0.0)
    end
    ω_grid, ω_w = gauleg(ωa, ωb, omega_nodes)
    phases = [_propagator_phase(meson, ω, q, qp, tp, K_coeffs; eta=eta) for ω in ω_grid]
    delta = _unwrap_phases(phases)
    F_current = delta
    F_gbu = _gbu_phase_function.(delta)

    current = 0.0
    gbu = 0.0
    @inbounds for iω in eachindex(ω_grid, ω_w, F_current, F_gbu)
        gω = bose_distribution(Float64(ω_grid[iω]), 0.0, Float64(tp.T))
        weight = ω_w[iω] * gω * (1.0 + gω) / (2.0 * π * Float64(tp.T))
        current += weight * F_current[iω]
        gbu += weight * F_gbu[iω]
    end
    return (current=current, gbu=gbu)
end

function _threshold_split_components(meson::Symbol, threshold_ref::Float64, qp, tp, K_coeffs;
                                     degeneracy::Int, qmax::Float64, q_nodes::Int,
                                     omega_min::Float64, omega_max::Float64, omega_nodes::Int,
                                     eta::Float64)
    q_grid, q_w = gauleg(0.0, qmax, q_nodes)

    current_below_sum = 0.0
    current_above_sum = 0.0
    gbu_below_sum = 0.0
    gbu_above_sum = 0.0

    @inbounds for iq in eachindex(q_grid, q_w)
        q = q_grid[iq]
        ωthr = hypot(q, threshold_ref)

        below_hi = min(ωthr, omega_max)
        above_lo = max(ωthr, omega_min)

        below = _omega_phase_integral(
            meson, q, omega_min, below_hi, qp, tp, K_coeffs;
            omega_nodes=omega_nodes, eta=eta,
        )
        above = _omega_phase_integral(
            meson, q, above_lo, omega_max, qp, tp, K_coeffs;
            omega_nodes=omega_nodes, eta=eta,
        )

        q_pref = q^2 / (2.0 * π^2)
        current_below_sum += q_w[iq] * q_pref * below.current
        current_above_sum += q_w[iq] * q_pref * above.current
        gbu_below_sum += q_w[iq] * q_pref * below.gbu
        gbu_above_sum += q_w[iq] * q_pref * above.gbu
    end

    pref = Float64(degeneracy)
    current_total = pref * (current_below_sum + current_above_sum)
    gbu_total = pref * (gbu_below_sum + gbu_above_sum)

    return (
        density_current_below=pref * current_below_sum,
        density_current_above=pref * current_above_sum,
        density_current_total=current_total,
        density_gbu_below=pref * gbu_below_sum,
        density_gbu_above=pref * gbu_above_sum,
        density_gbu_total=gbu_total,
        current_above_fraction=iszero(current_total) ? NaN : (pref * current_above_sum) / current_total,
        gbu_above_fraction=iszero(gbu_total) ? NaN : (pref * gbu_above_sum) / gbu_total,
    )
end

function main()
    output = isempty(ARGS) ? DEFAULT_OUTPUT : abspath(ARGS[1])
    mkpath(dirname(output))

    Ts_MeV = _selected_temperatures()
    qmax = parse(Float64, get(ENV, "MESON_E5_QMAX", string(DEFAULT_QMAX)))
    q_nodes = parse(Int, get(ENV, "MESON_E5_Q_NODES", string(DEFAULT_Q_NODES)))
    omega_min = parse(Float64, get(ENV, "MESON_E5_OMEGA_MIN", string(DEFAULT_OMEGA_MIN)))
    omega_max = parse(Float64, get(ENV, "MESON_E5_OMEGA_MAX", string(DEFAULT_OMEGA_MAX)))
    omega_nodes = parse(Int, get(ENV, "MESON_E5_OMEGA_NODES", string(DEFAULT_OMEGA_NODES)))
    eta = parse(Float64, get(ENV, "MESON_E5_ETA", string(DEFAULT_ETA)))

    continuation_state = nothing
    rows = NamedTuple[]

    for T_MeV in Ts_MeV
        T_fm = T_MeV / ħc_MeV_fm
        meson_point = solve_gap_and_meson_point(
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
        continuation_state = meson_point.continuation_state

        tp = normalize_thermo_params(meson_point.thermo_params)
        qp = ensure_quark_params_has_A(normalize_quark_params(meson_point.quark_params), tp)
        K_coeffs = _build_k_coeffs(qp)

        for meson in (:pi, :K)
            mr = meson_point.meson_results[meson]
            split = _threshold_split_components(
                meson,
                Float64(mr.threshold),
                qp,
                tp,
                K_coeffs;
                degeneracy=meson_degeneracy(meson),
                qmax=qmax,
                q_nodes=q_nodes,
                omega_min=omega_min,
                omega_max=omega_max,
                omega_nodes=omega_nodes,
                eta=eta,
            )

            push!(rows, (
                T_MeV=T_MeV,
                meson=String(meson),
                mass_ref=Float64(mr.mass),
                gamma_ref=Float64(mr.gamma),
                threshold_ref=Float64(mr.threshold),
                density_current_below=split.density_current_below,
                density_current_above=split.density_current_above,
                density_current_total=split.density_current_total,
                density_gbu_below=split.density_gbu_below,
                density_gbu_above=split.density_gbu_above,
                density_gbu_total=split.density_gbu_total,
                current_above_fraction=split.current_above_fraction,
                gbu_above_fraction=split.gbu_above_fraction,
                qmax=qmax,
                q_nodes=q_nodes,
                omega_min=omega_min,
                omega_max=omega_max,
                omega_nodes=omega_nodes,
                eta=eta,
            ))
        end
    end

    open(output, "w") do io
        println(io, "T_MeV,meson,mass_ref,gamma_ref,threshold_ref,density_current_below,density_current_above,density_current_total,density_gbu_below,density_gbu_above,density_gbu_total,current_above_fraction,gbu_above_fraction,qmax,q_nodes,omega_min,omega_max,omega_nodes,eta")
        for r in rows
            println(io, join((_fmt(getproperty(r, k)) for k in (
                :T_MeV, :meson, :mass_ref, :gamma_ref, :threshold_ref,
                :density_current_below, :density_current_above, :density_current_total,
                :density_gbu_below, :density_gbu_above, :density_gbu_total,
                :current_above_fraction, :gbu_above_fraction,
                :qmax, :q_nodes, :omega_min, :omega_max, :omega_nodes, :eta
            )), ','))
        end
    end

    println("Wrote Phase E5 threshold-split CSV: ", output)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
