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
- 有限窗口比较显式包含 `[g(omega)F(omega)]_{omega_min}^{omega_max}`
- 端点和 GL 内点在同一个序列上 unwrap；本门禁只认证 smooth/unwrapped branch
- 旧 `density_*_byparts` / `*_diff` 列保留为 bulk-only 历史映射，closure 列才是恒等式门禁
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "workflow_apps", "WorkflowParamAdapters.jl"))
include(joinpath(@__DIR__, "bu_kernel_gate_utils.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .Models: solve_gap_and_meson_point
using .WorkflowParamAdapters: normalize_quark_params, normalize_thermo_params
using Main.GaussLegendre: gauleg
using Main.MesonDensity: bose_distribution, meson_degeneracy
using Main.MesonMass: ensure_quark_params_has_A
using Main.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using Main.PolarizationAniso: polarization_with_width
using Main.MesonPropagator: meson_propagator_simple
using .BUKernelGateUtils: finite_window_bu_identity, nonuniform_three_point_derivative

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

function _channel_consistency(meson::Symbol, qp, tp, K_coeffs;
                              degeneracy::Int, qmax::Float64, q_nodes::Int,
                              omega_min::Float64, omega_max::Float64, omega_nodes::Int,
                              eta::Float64)
    q_grid, q_w = gauleg(0.0, qmax, q_nodes)
    omega_grid, omega_w = gauleg(omega_min, omega_max, omega_nodes)

    current_bulk_sum = 0.0
    current_boundary_lower_sum = 0.0
    current_boundary_upper_sum = 0.0
    current_derivative_sum = 0.0
    current_weighted_phase_min_sum = 0.0
    current_weighted_phase_max_sum = 0.0
    gbu_bulk_sum = 0.0
    gbu_boundary_lower_sum = 0.0
    gbu_boundary_upper_sum = 0.0
    gbu_derivative_sum = 0.0
    gbu_weighted_phase_min_sum = 0.0
    gbu_weighted_phase_max_sum = 0.0
    max_q_current_closure_abs = 0.0
    max_q_gbu_closure_abs = 0.0
    bose_omega_min = NaN
    bose_omega_max = NaN

    @inbounds for iq in eachindex(q_grid, q_w)
        q = q_grid[iq]
        omega_all = vcat(omega_min, collect(Float64, omega_grid), omega_max)
        phases_all = [
            _propagator_phase(meson, ω, q, qp, tp, K_coeffs; eta=eta)
            for ω in omega_all
        ]
        delta_all = _unwrap_phases(phases_all)
        ddelta_all = nonuniform_three_point_derivative(omega_all, delta_all)
        F_current_all = delta_all
        F_gbu_all = _gbu_phase_function.(delta_all)
        dF_gbu_all = _gbu_derivative_weight.(delta_all) .* ddelta_all

        interior = 2:(length(omega_all) - 1)
        current_gate = finite_window_bu_identity(
            omega_grid, omega_w, view(F_current_all, interior), view(ddelta_all, interior);
            T=Float64(tp.T), mu=0.0,
            omega_min=omega_min, omega_max=omega_max,
            F_min=first(F_current_all), F_max=last(F_current_all),
        )
        gbu_gate = finite_window_bu_identity(
            omega_grid, omega_w, view(F_gbu_all, interior), view(dF_gbu_all, interior);
            T=Float64(tp.T), mu=0.0,
            omega_min=omega_min, omega_max=omega_max,
            F_min=first(F_gbu_all), F_max=last(F_gbu_all),
        )

        q_pref = q^2 / (2.0 * π^2)
        q_weight = q_w[iq] * q_pref
        current_bulk_sum += q_weight * current_gate.byparts_bulk
        current_boundary_lower_sum += q_weight * current_gate.boundary_lower
        current_boundary_upper_sum += q_weight * current_gate.boundary_upper
        current_derivative_sum += q_weight * current_gate.derivative
        current_weighted_phase_min_sum += q_weight * current_gate.weighted_phase_min
        current_weighted_phase_max_sum += q_weight * current_gate.weighted_phase_max
        gbu_bulk_sum += q_weight * gbu_gate.byparts_bulk
        gbu_boundary_lower_sum += q_weight * gbu_gate.boundary_lower
        gbu_boundary_upper_sum += q_weight * gbu_gate.boundary_upper
        gbu_derivative_sum += q_weight * gbu_gate.derivative
        gbu_weighted_phase_min_sum += q_weight * gbu_gate.weighted_phase_min
        gbu_weighted_phase_max_sum += q_weight * gbu_gate.weighted_phase_max
        max_q_current_closure_abs = max(max_q_current_closure_abs, current_gate.closure_abs)
        max_q_gbu_closure_abs = max(max_q_gbu_closure_abs, gbu_gate.closure_abs)
        bose_omega_min = current_gate.bose_min
        bose_omega_max = current_gate.bose_max
    end

    pref = Float64(degeneracy)
    current_bulk = pref * current_bulk_sum
    current_boundary_lower = pref * current_boundary_lower_sum
    current_boundary_upper = pref * current_boundary_upper_sum
    current_boundary = current_boundary_lower + current_boundary_upper
    current_total = current_bulk + current_boundary
    current_derivative = pref * current_derivative_sum
    gbu_bulk = pref * gbu_bulk_sum
    gbu_boundary_lower = pref * gbu_boundary_lower_sum
    gbu_boundary_upper = pref * gbu_boundary_upper_sum
    gbu_boundary = gbu_boundary_lower + gbu_boundary_upper
    gbu_total = gbu_bulk + gbu_boundary
    gbu_derivative = pref * gbu_derivative_sum
    current_closure = current_derivative - current_total
    gbu_closure = gbu_derivative - gbu_total

    return (
        # Backward-compatible aliases: these remain the old bulk-only quantity.
        density_current_byparts=current_bulk,
        density_current_derivative=current_derivative,
        density_gbu_byparts=gbu_bulk,
        density_gbu_derivative=gbu_derivative,
        current_abs_diff=abs(current_derivative - current_bulk),
        current_rel_diff=abs(current_derivative - current_bulk) / max(abs(current_bulk), 1e-12),
        gbu_abs_diff=abs(gbu_derivative - gbu_bulk),
        gbu_rel_diff=abs(gbu_derivative - gbu_bulk) / max(abs(gbu_bulk), 1e-12),
        density_current_byparts_bulk=current_bulk,
        density_current_boundary_lower=current_boundary_lower,
        density_current_boundary_upper=current_boundary_upper,
        density_current_boundary=current_boundary,
        density_current_byparts_total=current_total,
        current_closure_abs=abs(current_closure),
        current_closure_rel=abs(current_closure) / max(abs(current_derivative), abs(current_total), 1e-12),
        density_gbu_byparts_bulk=gbu_bulk,
        density_gbu_boundary_lower=gbu_boundary_lower,
        density_gbu_boundary_upper=gbu_boundary_upper,
        density_gbu_boundary=gbu_boundary,
        density_gbu_byparts_total=gbu_total,
        gbu_closure_abs=abs(gbu_closure),
        gbu_closure_rel=abs(gbu_closure) / max(abs(gbu_derivative), abs(gbu_total), 1e-12),
        current_qweighted_phase_omega_min=pref * current_weighted_phase_min_sum,
        current_qweighted_phase_omega_max=pref * current_weighted_phase_max_sum,
        gbu_qweighted_phase_omega_min=pref * gbu_weighted_phase_min_sum,
        gbu_qweighted_phase_omega_max=pref * gbu_weighted_phase_max_sum,
        bose_omega_min=bose_omega_min,
        bose_omega_max=bose_omega_max,
        max_q_current_closure_abs=max_q_current_closure_abs,
        max_q_gbu_closure_abs=max_q_gbu_closure_abs,
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
                density_current_byparts_bulk=cc.density_current_byparts_bulk,
                density_current_boundary_lower=cc.density_current_boundary_lower,
                density_current_boundary_upper=cc.density_current_boundary_upper,
                density_current_boundary=cc.density_current_boundary,
                density_current_byparts_total=cc.density_current_byparts_total,
                current_closure_abs=cc.current_closure_abs,
                current_closure_rel=cc.current_closure_rel,
                density_gbu_byparts_bulk=cc.density_gbu_byparts_bulk,
                density_gbu_boundary_lower=cc.density_gbu_boundary_lower,
                density_gbu_boundary_upper=cc.density_gbu_boundary_upper,
                density_gbu_boundary=cc.density_gbu_boundary,
                density_gbu_byparts_total=cc.density_gbu_byparts_total,
                gbu_closure_abs=cc.gbu_closure_abs,
                gbu_closure_rel=cc.gbu_closure_rel,
                current_qweighted_phase_omega_min=cc.current_qweighted_phase_omega_min,
                current_qweighted_phase_omega_max=cc.current_qweighted_phase_omega_max,
                gbu_qweighted_phase_omega_min=cc.gbu_qweighted_phase_omega_min,
                gbu_qweighted_phase_omega_max=cc.gbu_qweighted_phase_omega_max,
                bose_omega_min=cc.bose_omega_min,
                bose_omega_max=cc.bose_omega_max,
                max_q_current_closure_abs=cc.max_q_current_closure_abs,
                max_q_gbu_closure_abs=cc.max_q_gbu_closure_abs,
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
        println(io, "T_MeV,meson,mass_ref,gamma_ref,threshold_ref,density_current_byparts,density_current_derivative,density_gbu_byparts,density_gbu_derivative,current_abs_diff,current_rel_diff,gbu_abs_diff,gbu_rel_diff,density_current_byparts_bulk,density_current_boundary_lower,density_current_boundary_upper,density_current_boundary,density_current_byparts_total,current_closure_abs,current_closure_rel,density_gbu_byparts_bulk,density_gbu_boundary_lower,density_gbu_boundary_upper,density_gbu_boundary,density_gbu_byparts_total,gbu_closure_abs,gbu_closure_rel,current_qweighted_phase_omega_min,current_qweighted_phase_omega_max,gbu_qweighted_phase_omega_min,gbu_qweighted_phase_omega_max,bose_omega_min,bose_omega_max,max_q_current_closure_abs,max_q_gbu_closure_abs,qmax,q_nodes,omega_min,omega_max,omega_nodes,eta")
        for r in rows
            println(io, join((_fmt(getproperty(r, k)) for k in (
                :T_MeV, :meson, :mass_ref, :gamma_ref, :threshold_ref,
                :density_current_byparts, :density_current_derivative,
                :density_gbu_byparts, :density_gbu_derivative,
                :current_abs_diff, :current_rel_diff,
                :gbu_abs_diff, :gbu_rel_diff,
                :density_current_byparts_bulk,
                :density_current_boundary_lower, :density_current_boundary_upper,
                :density_current_boundary, :density_current_byparts_total,
                :current_closure_abs, :current_closure_rel,
                :density_gbu_byparts_bulk,
                :density_gbu_boundary_lower, :density_gbu_boundary_upper,
                :density_gbu_boundary, :density_gbu_byparts_total,
                :gbu_closure_abs, :gbu_closure_rel,
                :current_qweighted_phase_omega_min, :current_qweighted_phase_omega_max,
                :gbu_qweighted_phase_omega_min, :gbu_qweighted_phase_omega_max,
                :bose_omega_min, :bose_omega_max,
                :max_q_current_closure_abs, :max_q_gbu_closure_abs,
                :qmax, :q_nodes, :omega_min, :omega_max, :omega_nodes, :eta
            )), ','))
        end
    end

    println("Wrote Phase E5 derivative-vs-byparts CSV: ", output)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
