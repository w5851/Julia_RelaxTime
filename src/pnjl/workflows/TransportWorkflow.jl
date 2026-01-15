module TransportWorkflow

"""
将"各向异性 PNJL 能隙方程求解"与"RTA 输运系数计算"串联起来，并可选在内部计算弛豫时间 τ。

工作流：
1) 调用 `PNJL.solve(FixedMu(), T_fm, μ_fm)` 在给定 `(T, μ, ξ)` 下求平衡（能隙/Polyakov 环）。
2) 从平衡解构造 `quark_params` 与 `thermo_params`。
3) （可选）计算数密度与平均散射率，从而得到 `tau`（见 `RelaxationTime.relaxation_times`）。
4) （可选）通过 `ThermoDerivatives.bulk_viscosity_coefficients` 生成体粘滞需要的等熵系数。
5) 调用 `TransportCoefficients.transport_coefficients` 返回 (η, ζ, σ)。

单位约定：
- `T_fm`, `mu_fm`, 质量/动量均为 fm⁻¹
- `tau` 为 fm

注意：仓库目前采用 include 组织方式，本模块会 include 相关依赖文件。
"""

# --- 依赖：常量/PNJL模块/弛豫时间/输运系数/一圈积分(A) ---
include("../../Constants_PNJL.jl")
include("../PNJL.jl")
include("../../relaxtime/RelaxationTime.jl")
include("../../relaxtime/TransportCoefficients.jl")
include("../../relaxtime/OneLoopIntegrals.jl")

using StaticArrays

using .PNJL: solve, FixedMu, cached_nodes, calculate_mass_vec, calculate_number_densities
using .PNJL: HADRON_SEED_5, DEFAULT_MOMENTUM_COUNT, DEFAULT_THETA_COUNT
using .PNJL.ThermoDerivatives: bulk_viscosity_coefficients
using .PNJL.Integrals: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .RelaxationTime: relaxation_times
using .TransportCoefficients: transport_coefficients
using .OneLoopIntegrals: A

export solve_gap_and_transport, build_equilibrium_params

"""将平衡求解结果转换成 (quark_params, thermo_params)。"""
function build_equilibrium_params(base, T_fm::Real, mu_fm::Real; xi::Real=0.0)
    Φ = Float64(base.x_state[4])
    Φbar = Float64(base.x_state[5])
    return (
        quark_params=(
            m=(u=NaN, d=NaN, s=NaN),
            μ=(u=Float64(mu_fm), d=Float64(mu_fm), s=Float64(mu_fm)),
        ),
        thermo_params=(T=Float64(T_fm), Φ=Φ, Φbar=Φbar, ξ=Float64(xi)),
    )
end

@inline function _densities_from_equilibrium(x_state, mu_vec, T_fm, thermal_nodes, xi)
    nd = calculate_number_densities(x_state, mu_vec, T_fm, thermal_nodes, xi)
    return (
        u=Float64(nd.quark[1]),
        d=Float64(nd.quark[2]),
        s=Float64(nd.quark[3]),
        ubar=Float64(nd.antiquark[1]),
        dbar=Float64(nd.antiquark[2]),
        sbar=Float64(nd.antiquark[3]),
    )
end

@inline function _A_from_equilibrium(T_fm::Real, quark_params::NamedTuple, thermo_params::NamedTuple)
    μu = quark_params.μ.u
    μs = quark_params.μ.s
    mu = quark_params.m.u
    ms = quark_params.m.s
    Φ = thermo_params.Φ
    Φbar = thermo_params.Φbar

    nodes = DEFAULT_MOMENTUM_NODES
    weights = DEFAULT_MOMENTUM_WEIGHTS

    A_u = A(Float64(mu), Float64(μu), Float64(T_fm), Φ, Φbar, nodes, weights)
    A_s = A(Float64(ms), Float64(μs), Float64(T_fm), Φ, Φbar, nodes, weights)
    return (u=A_u, d=A_u, s=A_s)
end

"""一次性完成：平衡求解 →（可选）τ 计算 →（可选）ζ 导数 → 输运系数。

返回值（NamedTuple）包含：
- `equilibrium`: 求解输出（pressure/energy/rho/.../x_state 等）
- `quark_params`, `thermo_params`
- `masses`: 三味有效质量
- `densities`: 用于 τ 计算的六种粒子/反粒子数密度
- `tau`, `tau_inv`: 弛豫时间及其倒数（若内部计算或调用方提供）
- `rates`: 平均散射率（若内部计算 τ 则返回，便于复用/诊断）
- `bulk_coeffs`: 若 `compute_bulk=true`，则给出体粘滞需要的导数组合
- `transport`: (eta, zeta, sigma)

关键词参数：
- `xi`: 各向异性参数 ξ
- `compute_tau`: 是否在内部计算 τ（需要 `K_coeffs`）
- `K_coeffs`: τ 所需的有效耦合系数（传给截面/平均散射率）
- `tau`: 若提供则直接用于输运计算，并仍会被原样回传（适合你自己在外部算 τ 或做对照）
- `compute_bulk`: 是否计算 ζ 需要的导数（会额外触发多次自动微分+求解，较慢）
- `p_num`, `t_num`: 热积分节点数（传给能隙求解/密度计算）
- `solver_kwargs`: 透传到 `solve`（例如 `iterations` 等）
- `tau_kwargs`: 透传到 `RelaxationTime.relaxation_times`（例如 `p_nodes/angle_nodes/phi_nodes/n_sigma_points/cs_caches/existing_rates` 等）
- `transport_kwargs`: 透传到 `transport_coefficients`（例如 `p_nodes`, `p_max` 等）
"""
function solve_gap_and_transport(
    T_fm::Real,
    mu_fm::Real;
    xi::Real=0.0,
    equilibrium::Union{Nothing,Any}=nothing,
    compute_tau::Bool=false,
    K_coeffs::Union{Nothing,NamedTuple}=nothing,
    tau::Union{Nothing,NamedTuple}=nothing,
    compute_bulk::Bool=true,
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    seed_state=HADRON_SEED_5,
    solver_kwargs::NamedTuple=(;),
    tau_kwargs::NamedTuple=(;),
    transport_kwargs::NamedTuple=(;)
)
    base = equilibrium === nothing ? begin
        # 兼容旧接口：允许直接传入一个状态向量作为初值。
        # 若调用方想使用更复杂的初值逻辑（MultiSeed/PhaseAwareContinuitySeed 等），应直接在外部调用 PNJL.solve。
        seed_strategy = if seed_state isa AbstractVector
            s5 = Float64.(seed_state[1:5])
            PNJL.DefaultSeed(s5, s5, :hadron)
        else
            PNJL.DefaultSeed(phase_hint=:auto)
        end

        solve(FixedMu(), T_fm, mu_fm;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            seed_strategy=seed_strategy,
            solver_kwargs...,
        )
    end : equilibrium

    # 有效质量（直接由平衡解得到）
    masses = base.masses

    params0 = build_equilibrium_params(base, T_fm, mu_fm; xi=xi)
    quark_params_basic = (
        m=(u=Float64(masses[1]), d=Float64(masses[2]), s=Float64(masses[3])),
        μ=params0.quark_params.μ,
    )
    thermo_params = params0.thermo_params

    # 密度：用于 τ 的 ω_i = Σ ρ_j \bar{w}_{ij}
    thermal_nodes = cached_nodes(p_num, t_num)
    densities = _densities_from_equilibrium(base.x_state, base.mu_vec, T_fm, thermal_nodes, Float64(xi))

    tau_inv = nothing
    rates = nothing

    if compute_tau
        K_coeffs === nothing && error("compute_tau=true requires K_coeffs")

        # 为截面/传播子准备 A 字段（TotalPropagator 会用到）
        A_vals = _A_from_equilibrium(T_fm, quark_params_basic, thermo_params)
        quark_params_full = (m=quark_params_basic.m, μ=quark_params_basic.μ, A=A_vals)

        tau_res = relaxation_times(
            quark_params_full,
            thermo_params,
            K_coeffs;
            densities=densities,
            tau_kwargs...,
        )
        tau = tau_res.tau
        tau_inv = tau_res.tau_inv
        rates = tau_res.rates
    else
        # 若不内部计算 τ，则要求调用方提供 τ
        tau === nothing && error("either provide tau=... or set compute_tau=true")
    end

    bulk_coeffs = nothing
    if compute_bulk
        # ζ 需要 ThermoDerivatives.bulk_viscosity_coefficients 的整套等熵系数（v_n_sq, dμB/dT|σ, dM/dT, dM/dμB, ...）
        # 该函数不依赖 seed_state，也不接受求解器的 iterations 等参数。
        bulk_coeffs = bulk_viscosity_coefficients(
            T_fm,
            mu_fm;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
        )
    end

    tr = transport_coefficients(
        quark_params_basic,
        thermo_params;
        tau=tau,
        bulk_coeffs=bulk_coeffs,
        transport_kwargs...,
    )

    return (
        equilibrium=base,
        quark_params=quark_params_basic,
        thermo_params=thermo_params,
        masses=masses,
        densities=densities,
        tau=tau,
        tau_inv=tau_inv,
        rates=rates,
        bulk_coeffs=bulk_coeffs,
        transport=tr,
    )
end

end # module
