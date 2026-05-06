"""
Phase E5: 更严格 BU 收口的最小候选对照试算。

目标：
- 固定与 Phase-E3 正式扫描一致的 meson workflow / continuation 口径
- 在同一组 `(T, q, omega)` 网格上，对比若干相位函数候选：
  1. 当前 Phase-E3: `F(delta) = delta`
  2. generalized BU 候选: `F(delta) = delta - 0.5 * sin(2delta)`
  3. 诊断性尾部平移: `delta -> delta - delta(omega_max, q)`
- 先判断数量级差异是否主要来自阈值以上宽能区的正相位面积

注意：
- 这是分析脚本，不直接改正式 helper / workflow 合同
- 尾部平移只作诊断，不自动视为最终物理口径
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "workflow_apps", "WorkflowParamAdapters.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .Models: solve_gap_and_meson_point, solve_meson_density_from_meson_point
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
    "meson_density_phase_e5_strict_candidates.csv",
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

@inline _gbu_phase_function(δ::Float64) = δ - 0.5 * sin(2.0 * δ)

function _omega_integral_for_phase(phase_values::Vector{Float64}, omega_grid, omega_w, T::Float64)
    total = 0.0
    @inbounds for iω in eachindex(omega_grid, omega_w, phase_values)
        gω = bose_distribution(Float64(omega_grid[iω]), 0.0, T)
        total += omega_w[iω] * gω * (1.0 + gω) * phase_values[iω]
    end
    return total / (2.0 * π)
end

function _candidate_density_components(meson::Symbol, qp, tp, K_coeffs;
                                       degeneracy::Int,
                                       qmax::Float64,
                                       q_nodes::Int,
                                       omega_min::Float64,
                                       omega_max::Float64,
                                       omega_nodes::Int,
                                       eta::Float64)
    q_grid, q_w = gauleg(0.0, qmax, q_nodes)
    omega_grid, omega_w = gauleg(omega_min, omega_max, omega_nodes)

    current_sum = 0.0
    gbu_sum = 0.0
    shifted_sum = 0.0
    gbu_shifted_sum = 0.0
    current_q0 = NaN
    gbu_q0 = NaN
    shifted_q0 = NaN
    gbu_shifted_q0 = NaN

    @inbounds for iq in eachindex(q_grid, q_w)
        q = q_grid[iq]
        phases = [_propagator_phase(meson, ω, q, qp, tp, K_coeffs; eta=eta) for ω in omega_grid]
        phase_unwrapped = _unwrap_phases(phases)
        tail_ref = phase_unwrapped[end]
        phase_shifted = phase_unwrapped .- tail_ref
        phase_gbu = _gbu_phase_function.(phase_unwrapped)
        phase_gbu_shifted = _gbu_phase_function.(phase_shifted)

        omega_current = _omega_integral_for_phase(phase_unwrapped, omega_grid, omega_w, Float64(tp.T))
        omega_gbu = _omega_integral_for_phase(phase_gbu, omega_grid, omega_w, Float64(tp.T))
        omega_shifted = _omega_integral_for_phase(phase_shifted, omega_grid, omega_w, Float64(tp.T))
        omega_gbu_shifted = _omega_integral_for_phase(phase_gbu_shifted, omega_grid, omega_w, Float64(tp.T))

        q_pref = q^2 / (2.0 * π^2)
        current_shell = q_pref * omega_current
        gbu_shell = q_pref * omega_gbu
        shifted_shell = q_pref * omega_shifted
        gbu_shifted_shell = q_pref * omega_gbu_shifted

        current_sum += q_w[iq] * current_shell
        gbu_sum += q_w[iq] * gbu_shell
        shifted_sum += q_w[iq] * shifted_shell
        gbu_shifted_sum += q_w[iq] * gbu_shifted_shell

        if iq == 1
            current_q0 = current_shell
            gbu_q0 = gbu_shell
            shifted_q0 = shifted_shell
            gbu_shifted_q0 = gbu_shifted_shell
        end
    end

    pref = Float64(degeneracy) / Float64(tp.T)
    return (
        density_current=pref * current_sum,
        density_gbu=pref * gbu_sum,
        density_tail_shifted=pref * shifted_sum,
        density_gbu_tail_shifted=pref * gbu_shifted_sum,
        q_integral_current=current_sum,
        q_integral_gbu=gbu_sum,
        q_integral_tail_shifted=shifted_sum,
        q_integral_gbu_tail_shifted=gbu_shifted_sum,
        q0_shell_current=current_q0,
        q0_shell_gbu=gbu_q0,
        q0_shell_tail_shifted=shifted_q0,
        q0_shell_gbu_tail_shifted=gbu_shifted_q0,
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

        stable = solve_meson_density_from_meson_point(meson_point; num_q_nodes=256)
        tp = normalize_thermo_params(meson_point.thermo_params)
        qp = ensure_quark_params_has_A(normalize_quark_params(meson_point.quark_params), tp)
        K_coeffs = _build_k_coeffs(qp)

        for meson in (:pi, :K)
            mr = meson_point.meson_results[meson]
            parts = _candidate_density_components(
                meson, qp, tp, K_coeffs;
                degeneracy=meson_degeneracy(meson),
                qmax=qmax,
                q_nodes=q_nodes,
                omega_min=omega_min,
                omega_max=omega_max,
                omega_nodes=omega_nodes,
                eta=eta,
            )

            stable_density = meson === :pi ? stable.n_pi : stable.n_K
            push!(rows, (
                T_MeV=T_MeV,
                meson=String(meson),
                mass_ref=Float64(mr.mass),
                gamma_ref=Float64(mr.gamma),
                threshold_ref=Float64(mr.threshold),
                stable_density=stable_density,
                density_current=parts.density_current,
                density_gbu=parts.density_gbu,
                density_tail_shifted=parts.density_tail_shifted,
                density_gbu_tail_shifted=parts.density_gbu_tail_shifted,
                current_over_stable=iszero(stable_density) ? NaN : parts.density_current / stable_density,
                gbu_over_stable=iszero(stable_density) ? NaN : parts.density_gbu / stable_density,
                tail_shifted_over_stable=iszero(stable_density) ? NaN : parts.density_tail_shifted / stable_density,
                gbu_tail_shifted_over_stable=iszero(stable_density) ? NaN : parts.density_gbu_tail_shifted / stable_density,
                q_integral_current=parts.q_integral_current,
                q_integral_gbu=parts.q_integral_gbu,
                q_integral_tail_shifted=parts.q_integral_tail_shifted,
                q_integral_gbu_tail_shifted=parts.q_integral_gbu_tail_shifted,
                q0_shell_current=parts.q0_shell_current,
                q0_shell_gbu=parts.q0_shell_gbu,
                q0_shell_tail_shifted=parts.q0_shell_tail_shifted,
                q0_shell_gbu_tail_shifted=parts.q0_shell_gbu_tail_shifted,
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
        println(io, "T_MeV,meson,mass_ref,gamma_ref,threshold_ref,stable_density,density_current,density_gbu,density_tail_shifted,density_gbu_tail_shifted,current_over_stable,gbu_over_stable,tail_shifted_over_stable,gbu_tail_shifted_over_stable,q_integral_current,q_integral_gbu,q_integral_tail_shifted,q_integral_gbu_tail_shifted,q0_shell_current,q0_shell_gbu,q0_shell_tail_shifted,q0_shell_gbu_tail_shifted,qmax,q_nodes,omega_min,omega_max,omega_nodes,eta")
        for r in rows
            println(io, join((_fmt(getproperty(r, k)) for k in (
                :T_MeV, :meson, :mass_ref, :gamma_ref, :threshold_ref, :stable_density,
                :density_current, :density_gbu, :density_tail_shifted, :density_gbu_tail_shifted,
                :current_over_stable, :gbu_over_stable, :tail_shifted_over_stable, :gbu_tail_shifted_over_stable,
                :q_integral_current, :q_integral_gbu, :q_integral_tail_shifted, :q_integral_gbu_tail_shifted,
                :q0_shell_current, :q0_shell_gbu, :q0_shell_tail_shifted, :q0_shell_gbu_tail_shifted,
                :qmax, :q_nodes, :omega_min, :omega_max, :omega_nodes, :eta
            )), ','))
        end
    end

    println("Wrote Phase E5 strict-candidate CSV: ", output)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
