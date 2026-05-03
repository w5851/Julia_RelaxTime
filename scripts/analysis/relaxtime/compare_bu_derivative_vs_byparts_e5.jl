"""
Phase E5: 对比 BU 导数形式与分部积分形式的一致性。

目标：
- 固定与 Phase-E3 / E5 一致的 meson workflow / continuation 口径
- 在同一组 `(T, q, omega)` 网格上，同时计算：
  1. current: `F(delta) = delta`
  2. generalized BU: `F(delta) = delta - 0.5 sin(2 delta)`
- 对比它们在
  - 导数形式
  - 分部积分形式
  下的离散一致性

说明：
- 当前仍是分析脚本，不改正式 helper / workflow 契约
- 由于 `omega_min > 0` 且 `omega_max < inf`，两种离散形式不要求点对点完全一致
- 该脚本主要用于判断：差异是“有限窗口/离散误差”，还是“当前相位函数口径本身不稳”
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
    "meson_density_phase_e5_derivative_vs_byparts.csv",
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
@inline _gbu_derivative_weight(δ::Float64) = 2.0 * sin(δ)^2

function _nonuniform_derivative(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    n == length(y) || throw(ArgumentError("x/y length mismatch"))
    n >= 2 || throw(ArgumentError("need at least 2 points for derivative"))
    dy = similar(y)
    if n == 2
        slope = (y[2] - y[1]) / (x[2] - x[1])
        dy[1] = slope
        dy[2] = slope
        return dy
    end
    dy[1] = (y[2] - y[1]) / (x[2] - x[1])
    for i in 2:(n - 1)
        dy[i] = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
    end
    dy[n] = (y[n] - y[n - 1]) / (x[n] - x[n - 1])
    return dy
end

function _omega_integral_byparts(F_values::Vector{Float64}, omega_grid, omega_w, T::Float64)
    total = 0.0
    @inbounds for iω in eachindex(omega_grid, omega_w, F_values)
        gω = bose_distribution(Float64(omega_grid[iω]), 0.0, T)
        total += omega_w[iω] * gω * (1.0 + gω) * F_values[iω]
    end
    total / (2.0 * π * T)
end

function _omega_integral_derivative(delta_values::Vector{Float64}, ddelta_values::Vector{Float64},
                                    omega_grid, omega_w, T::Float64; generalized::Bool=false)
    total = 0.0
    @inbounds for iω in eachindex(omega_grid, omega_w, delta_values, ddelta_values)
        gω = bose_distribution(Float64(omega_grid[iω]), 0.0, T)
        weight = generalized ? _gbu_derivative_weight(delta_values[iω]) : 1.0
        total += omega_w[iω] * gω * weight * ddelta_values[iω]
    end
    total / (2.0 * π)
end

function _channel_consistency(meson::Symbol, qp, tp, K_coeffs;
                              degeneracy::Int, qmax::Float64, q_nodes::Int,
                              omega_min::Float64, omega_max::Float64, omega_nodes::Int,
                              eta::Float64)
    q_grid, q_w = gauleg(0.0, qmax, q_nodes)
    omega_grid, omega_w = gauleg(omega_min, omega_max, omega_nodes)

    current_byparts_sum = 0.0
    current_derivative_sum = 0.0
    gbu_byparts_sum = 0.0
    gbu_derivative_sum = 0.0

    @inbounds for iq in eachindex(q_grid, q_w)
        q = q_grid[iq]
        phases = [_propagator_phase(meson, ω, q, qp, tp, K_coeffs; eta=eta) for ω in omega_grid]
        delta = _unwrap_phases(phases)
        ddelta = _nonuniform_derivative(collect(Float64, omega_grid), delta)
        F_current = delta
        F_gbu = _gbu_phase_function.(delta)

        omega_current_byparts = _omega_integral_byparts(F_current, omega_grid, omega_w, Float64(tp.T))
        omega_current_derivative = _omega_integral_derivative(delta, ddelta, omega_grid, omega_w, Float64(tp.T))
        omega_gbu_byparts = _omega_integral_byparts(F_gbu, omega_grid, omega_w, Float64(tp.T))
        omega_gbu_derivative = _omega_integral_derivative(delta, ddelta, omega_grid, omega_w, Float64(tp.T); generalized=true)

        q_pref = q^2 / (2.0 * π^2)
        current_byparts_sum += q_w[iq] * q_pref * omega_current_byparts
        current_derivative_sum += q_w[iq] * q_pref * omega_current_derivative
        gbu_byparts_sum += q_w[iq] * q_pref * omega_gbu_byparts
        gbu_derivative_sum += q_w[iq] * q_pref * omega_gbu_derivative
    end

    pref = Float64(degeneracy)
    current_byparts = pref * current_byparts_sum
    current_derivative = pref * current_derivative_sum
    gbu_byparts = pref * gbu_byparts_sum
    gbu_derivative = pref * gbu_derivative_sum

    return (
        density_current_byparts=current_byparts,
        density_current_derivative=current_derivative,
        density_gbu_byparts=gbu_byparts,
        density_gbu_derivative=gbu_derivative,
        current_abs_diff=abs(current_derivative - current_byparts),
        current_rel_diff=abs(current_derivative - current_byparts) / max(abs(current_byparts), 1e-12),
        gbu_abs_diff=abs(gbu_derivative - gbu_byparts),
        gbu_rel_diff=abs(gbu_derivative - gbu_byparts) / max(abs(gbu_byparts), 1e-12),
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
            cc = _channel_consistency(
                meson, qp, tp, K_coeffs;
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
                density_current_byparts=cc.density_current_byparts,
                density_current_derivative=cc.density_current_derivative,
                density_gbu_byparts=cc.density_gbu_byparts,
                density_gbu_derivative=cc.density_gbu_derivative,
                current_abs_diff=cc.current_abs_diff,
                current_rel_diff=cc.current_rel_diff,
                gbu_abs_diff=cc.gbu_abs_diff,
                gbu_rel_diff=cc.gbu_rel_diff,
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
        println(io, "T_MeV,meson,mass_ref,gamma_ref,threshold_ref,density_current_byparts,density_current_derivative,density_gbu_byparts,density_gbu_derivative,current_abs_diff,current_rel_diff,gbu_abs_diff,gbu_rel_diff,qmax,q_nodes,omega_min,omega_max,omega_nodes,eta")
        for r in rows
            println(io, join((_fmt(getproperty(r, k)) for k in (
                :T_MeV, :meson, :mass_ref, :gamma_ref, :threshold_ref,
                :density_current_byparts, :density_current_derivative,
                :density_gbu_byparts, :density_gbu_derivative,
                :current_abs_diff, :current_rel_diff,
                :gbu_abs_diff, :gbu_rel_diff,
                :qmax, :q_nodes, :omega_min, :omega_max, :omega_nodes, :eta
            )), ','))
        end
    end

    println("Wrote Phase E5 derivative-vs-byparts CSV: ", output)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
