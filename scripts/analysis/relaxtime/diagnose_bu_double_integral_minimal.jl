"""
Phase E3: 各向同性 (xi = 0) 下的最小 BU 双积分试探。

目标：
- 直接基于传播子相位 `delta_M(omega, q) = arg D_M`
- 采用分部积分后的 `delta_M` 本体形式
- 在 xi = 0 时把三维 q 积分降成径向一维积分
- 先诊断 q + omega 双积分是否稳定，不作为正式产出标准

当前固定范围：
- muB = 0
- xi = 0
- mesons = pi, K
- 少量温点
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "workflow_apps", "WorkflowParamAdapters.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .Models: solve_gap_and_meson_point
using .WorkflowParamAdapters: normalize_quark_params, normalize_thermo_params
using .GaussLegendre: gauleg
using Main.MesonMass: ensure_quark_params_has_A
using Main.MesonDensity: phase_shift_meson_density_summary

const DEFAULT_OUTPUT = joinpath(
    PROJECT_ROOT,
    "data",
    "outputs",
    "results",
    "relaxtime",
    "analysis",
    "bu_phase_shift_minimal",
    "bu_double_integral_minimal.csv",
)

const DEFAULT_T_VALUES_MEV = [208.0, 210.0, 212.0]
const DEFAULT_Q_MIN = 0.0
const DEFAULT_Q_MAX = 12.0
const DEFAULT_OMEGA_MIN = 0.05
const DEFAULT_OMEGA_MAX = 10.0
const DEFAULT_Q_NODES = 48
const DEFAULT_OMEGA_NODES = 48
const DEFAULT_ETA = 1e-6

@inline function _fmt(x)
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
    raw = strip(get(ENV, "PHASE_E_T_VALUES", ""))
    isempty(raw) && return DEFAULT_T_VALUES_MEV
    return _parse_float_list(raw)
end

function _q_nodes_weights()
    qmin = parse(Float64, get(ENV, "PHASE_E_Q_MIN", string(DEFAULT_Q_MIN)))
    qmax = parse(Float64, get(ENV, "PHASE_E_Q_MAX", string(DEFAULT_Q_MAX)))
    qnodes = parse(Int, get(ENV, "PHASE_E_Q_NODES", string(DEFAULT_Q_NODES)))
    qnodes > 1 || throw(ArgumentError("PHASE_E_Q_NODES must exceed 1"))
    qmax > qmin || throw(ArgumentError("PHASE_E_Q_MAX must exceed PHASE_E_Q_MIN"))
    return gauleg(qmin, qmax, qnodes)
end

function _omega_nodes_weights()
    ωmin = parse(Float64, get(ENV, "PHASE_E_OMEGA_MIN", string(DEFAULT_OMEGA_MIN)))
    ωmax = parse(Float64, get(ENV, "PHASE_E_OMEGA_MAX", string(DEFAULT_OMEGA_MAX)))
    ωnodes = parse(Int, get(ENV, "PHASE_E_OMEGA_NODES", string(DEFAULT_OMEGA_NODES)))
    ωnodes > 1 || throw(ArgumentError("PHASE_E_OMEGA_NODES must exceed 1"))
    ωmax > ωmin || throw(ArgumentError("PHASE_E_OMEGA_MAX must exceed PHASE_E_OMEGA_MIN"))
    return gauleg(ωmin, ωmax, ωnodes)
end

function q_nodes_weights_from_env()
    return _q_nodes_weights()
end

function omega_nodes_weights_from_env()
    return _omega_nodes_weights()
end

function solve_equilibrium_meson_point(T_MeV::Float64; continuation_state=nothing)
    T_fm = T_MeV / ħc_MeV_fm
    return solve_gap_and_meson_point(
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
end

function compute_double_integral_rows(
    equilibrium_result;
    q_nodes::AbstractVector{<:Real},
    q_weights::AbstractVector{<:Real},
    omega_nodes::AbstractVector{<:Real},
    omega_weights::AbstractVector{<:Real},
    qmax_request::Float64,
    omega_max_request::Float64,
    eta::Float64=DEFAULT_ETA,
)
    tp = normalize_thermo_params(equilibrium_result.thermo_params)
    qp = ensure_quark_params_has_A(normalize_quark_params(equilibrium_result.quark_params), tp)
    summary = phase_shift_meson_density_summary(
        qp,
        tp;
        qmax=qmax_request,
        q_nodes=length(q_nodes),
        omega_min=Float64(omega_nodes[1]) > 0.0 ? parse(Float64, get(ENV, "PHASE_E_OMEGA_MIN", string(DEFAULT_OMEGA_MIN))) : DEFAULT_OMEGA_MIN,
        omega_max=omega_max_request,
        omega_nodes=length(omega_nodes),
        eta=eta,
    )
    T_MeV = Float64(tp.T) * ħc_MeV_fm
    rows = [
        (
            T_MeV=T_MeV,
            meson=:pi,
            qmax_request=qmax_request,
            omega_max_request=omega_max_request,
            q_last_node=Float64(q_nodes[end]),
            omega_last_node=Float64(omega_nodes[end]),
            q_nodes=length(q_nodes),
            omega_nodes=length(omega_nodes),
            q_integral_estimate=Float64(summary.pi_density.q_integral_estimate),
            omega_shell_at_qmax=Float64(summary.pi_density.omega_shell_at_qmax),
            double_integral_estimate=Float64(summary.n_pi),
        ),
        (
            T_MeV=T_MeV,
            meson=:K,
            qmax_request=qmax_request,
            omega_max_request=omega_max_request,
            q_last_node=Float64(q_nodes[end]),
            omega_last_node=Float64(omega_nodes[end]),
            q_nodes=length(q_nodes),
            omega_nodes=length(omega_nodes),
            q_integral_estimate=Float64(summary.k_density.q_integral_estimate),
            omega_shell_at_qmax=Float64(summary.k_density.omega_shell_at_qmax),
            double_integral_estimate=Float64(summary.n_K),
        ),
    ]
    return rows
end

function main()
    env_output = strip(get(ENV, "PHASE_E_OUTPUT", ""))
    output = if !isempty(ARGS)
        abspath(ARGS[1])
    elseif !isempty(env_output)
        abspath(env_output)
    else
        DEFAULT_OUTPUT
    end
    mkpath(dirname(output))

    Ts_MeV = _selected_temperatures()
    q_nodes, q_weights = _q_nodes_weights()
    omega_nodes, omega_weights = _omega_nodes_weights()
    qmax_request = parse(Float64, get(ENV, "PHASE_E_Q_MAX", string(DEFAULT_Q_MAX)))
    omega_max_request = parse(Float64, get(ENV, "PHASE_E_OMEGA_MAX", string(DEFAULT_OMEGA_MAX)))
    eta = parse(Float64, get(ENV, "PHASE_E_ETA", string(DEFAULT_ETA)))

    continuation_state = nothing

    open(output, "w") do io
        println(io, "T_MeV,meson,qmax_request,omega_max_request,q_last_node,omega_last_node,q_nodes,omega_nodes,q_integral_estimate,omega_shell_at_qmax,double_integral_estimate")

        for T_MeV in Ts_MeV
            res = solve_equilibrium_meson_point(T_MeV; continuation_state=continuation_state)
            continuation_state = res.continuation_state
            rows = compute_double_integral_rows(
                res;
                q_nodes=q_nodes,
                q_weights=q_weights,
                omega_nodes=omega_nodes,
                omega_weights=omega_weights,
                qmax_request=qmax_request,
                omega_max_request=omega_max_request,
                eta=eta,
            )
            for row in rows
                println(io, join((
                    _fmt(row.T_MeV),
                    String(row.meson),
                    _fmt(row.qmax_request),
                    _fmt(row.omega_max_request),
                    _fmt(row.q_last_node),
                    _fmt(row.omega_last_node),
                    _fmt(row.q_nodes),
                    _fmt(row.omega_nodes),
                    _fmt(row.q_integral_estimate),
                    _fmt(row.omega_shell_at_qmax),
                    _fmt(row.double_integral_estimate),
                ), ','))
            end
        end
    end

    println("Wrote double-integral diagnostic CSV: ", output)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
