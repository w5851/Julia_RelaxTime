module MesonDensityWorkflow

"""
将现有 meson workflow 返回值后处理为介子数密度。

职责边界：
1. 不重写平衡态求解；
2. 不重写介子质量 / Mott 求解；
3. 只消费 `solve_gap_and_meson_point` 的返回值，生成 `π/K` 稳定粒子极限数密度与 `K/π` 比值。

这保证介子数密度主线继续以 meson workflow 为唯一计算入口，而不是在脚本层平行拼接流程。
"""

const _RELAX_TIME_MODULE_PATH = normpath(joinpath(@__DIR__, "..", "..", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAX_TIME_MODULE_PATH)
end

using ..MesonMassWorkflow: solve_gap_and_meson_point
using Main.MesonDensity: DEFAULT_MESON_DENSITY_Q_NODES,
                         DEFAULT_PHASE_SHIFT_Q_MAX,
                         DEFAULT_PHASE_SHIFT_Q_NODES,
                         DEFAULT_PHASE_SHIFT_OMEGA_MIN,
                         DEFAULT_PHASE_SHIFT_OMEGA_MAX,
                         DEFAULT_PHASE_SHIFT_OMEGA_NODES,
                         meson_degeneracy,
                         phase_shift_point_diagnostic,
                         phase_shift_meson_density_derivative_reference_summary,
                         phase_shift_meson_density_summary,
                         strict_bw_meson_density_summary,
                         strict_bw_qpole_density_summary,
                         stable_meson_number_density
using ..WorkflowParamAdapters: normalize_quark_params, normalize_thermo_params
using Main.MesonMass: ensure_quark_params_has_A

export solve_meson_density_from_meson_point
export solve_gap_and_meson_density_point
export solve_strict_bw_meson_density_from_meson_point
export solve_gap_and_strict_bw_meson_density_point
export solve_phase_shift_meson_density_from_meson_point
export solve_gap_and_phase_shift_meson_density_point
export solve_phase_shift_point_diagnostic_from_meson_point
export solve_phase_shift_derivative_reference_from_meson_point

@inline function _require_result_field(result, field::Symbol)
    hasproperty(result, field) || throw(ArgumentError("meson workflow result missing required field: $(field)"))
    return getproperty(result, field)
end

@inline function _require_finite_mass(meson_results, meson::Symbol)::Float64
    haskey(meson_results, meson) || throw(ArgumentError("meson_results missing required channel $(meson)"))
    entry = meson_results[meson]
    hasproperty(entry, :mass) || throw(ArgumentError("meson_results[$(meson)] missing mass field"))
    mass = Float64(entry.mass)
    isfinite(mass) || throw(ArgumentError("meson_results[$(meson)].mass must be finite, got $(mass)"))
    mass >= 0.0 || throw(ArgumentError(
        "stable meson-density workflow requires nonnegative $(meson) mass, got $(mass); " *
        "this point is outside the current stable-particle validity window and should be treated by BW/BU extensions instead"
    ))
    return mass
end

@inline function _require_finite_gamma(meson_results, meson::Symbol)::Float64
    haskey(meson_results, meson) || throw(ArgumentError("meson_results missing required channel $(meson)"))
    entry = meson_results[meson]
    hasproperty(entry, :gamma) || throw(ArgumentError("meson_results[$(meson)] missing gamma field"))
    gamma = Float64(entry.gamma)
    isfinite(gamma) || throw(ArgumentError("meson_results[$(meson)].gamma must be finite, got $(gamma)"))
    gamma >= 0.0 || throw(ArgumentError("meson_results[$(meson)].gamma must be nonnegative, got $(gamma)"))
    return gamma
end

"""
    solve_meson_density_from_meson_point(meson_point; kwargs...) -> NamedTuple

后处理入口：消费 `solve_gap_and_meson_point` 的返回值，计算稳定粒子极限
`n_pi(T)`、`n_K(T)` 与 `K/π`。

约束：
- 当前默认只读取 `meson_results[:pi]` 与 `meson_results[:K]`
- 温度来自 `meson_point.thermo_params.T`
- 不在此层重复求平衡态或介子极点
"""
function solve_meson_density_from_meson_point(
    meson_point;
    pi_channel::Symbol=:pi,
    k_channel::Symbol=:K,
    μ_pi::Real=0.0,
    μ_K::Real=0.0,
    d_pi::Integer=meson_degeneracy(:pi),
    d_K::Integer=meson_degeneracy(:K),
    qmax_pi::Union{Nothing,Real}=nothing,
    qmax_K::Union{Nothing,Real}=nothing,
    num_q_nodes::Int=DEFAULT_MESON_DENSITY_Q_NODES,
)
    thermo_params = _require_result_field(meson_point, :thermo_params)
    meson_results = _require_result_field(meson_point, :meson_results)

    T_fm = Float64(thermo_params.T)
    xi = Float64(getproperty(thermo_params, :ξ))
    m_pi = _require_finite_mass(meson_results, pi_channel)
    m_K = _require_finite_mass(meson_results, k_channel)

    qmax_pi_fm = qmax_pi === nothing ? nothing : Float64(qmax_pi)
    qmax_K_fm = qmax_K === nothing ? nothing : Float64(qmax_K)
    μ_pi_fm = Float64(μ_pi)
    μ_K_fm = Float64(μ_K)

    n_pi = stable_meson_number_density(
        m_pi,
        T_fm;
        μ=μ_pi_fm,
        degeneracy=Int(d_pi),
        qmax=qmax_pi_fm,
        num_q_nodes=num_q_nodes,
    )
    n_K = stable_meson_number_density(
        m_K,
        T_fm;
        μ=μ_K_fm,
        degeneracy=Int(d_K),
        qmax=qmax_K_fm,
        num_q_nodes=num_q_nodes,
    )

    return (
        T_fm=T_fm,
        xi=xi,
        pi_channel=pi_channel,
        k_channel=k_channel,
        m_pi=m_pi,
        m_K=m_K,
        μ_pi=μ_pi_fm,
        μ_K=μ_K_fm,
        d_pi=Int(d_pi),
        d_K=Int(d_K),
        n_pi=n_pi,
        n_K=n_K,
        kpi_ratio=iszero(n_pi) ? NaN : n_K / n_pi,
        num_q_nodes=num_q_nodes,
        qmax_pi=qmax_pi_fm,
        qmax_K=qmax_K_fm,
    )
end

"""
    solve_gap_and_meson_density_point(T_fm, mu_fm; density_kwargs=(;), kwargs...) -> NamedTuple

完整入口：先走 meson workflow，再在其返回值上追加 `meson_density` 后处理结果。

设计目标：
- 让脚本层只调用一个统一 workflow 入口；
- 避免脚本重新组装平衡态、介子质量和数密度链路。
"""
function solve_gap_and_meson_density_point(
    T_fm::Real,
    mu_fm::Real;
    density_kwargs::NamedTuple=(;),
    kwargs...,
)
    meson_point = solve_gap_and_meson_point(T_fm, mu_fm; kwargs...)
    density = solve_meson_density_from_meson_point(meson_point; density_kwargs...)
    return merge(meson_point, (meson_density=density,))
end

"""
    solve_strict_bw_meson_density_from_meson_point(meson_point; kwargs...) -> NamedTuple

消费 meson workflow 返回值，执行 Stage-1 reduced strict BW 数密度后处理。

当前约束：
- 只消费 workflow 已返回的 `mass` / `gamma`
- 不在此层重求 `q` 依赖复极点
"""
function solve_strict_bw_meson_density_from_meson_point(
    meson_point;
    pi_channel::Symbol=:pi,
    k_channel::Symbol=:K,
    stage::Symbol=:stage1_reduced,
    μ_pi::Real=0.0,
    μ_K::Real=0.0,
    d_pi::Integer=meson_degeneracy(:pi),
    d_K::Integer=meson_degeneracy(:K),
    qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    omega_min::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MIN,
    omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    gamma_zero_tol::Float64=1e-12,
    solver_iterations::Int=20,
    pole_residual_norm_max::Float64=1e-6,
    pole_require_converged::Bool=true,
)
    thermo_params = _require_result_field(meson_point, :thermo_params)
    meson_results = _require_result_field(meson_point, :meson_results)

    T_fm = Float64(thermo_params.T)
    xi = Float64(getproperty(thermo_params, :ξ))
    m_pi = _require_finite_mass(meson_results, pi_channel)
    m_K = _require_finite_mass(meson_results, k_channel)
    gamma_pi = _require_finite_gamma(meson_results, pi_channel)
    gamma_K = _require_finite_gamma(meson_results, k_channel)

    stage_sym = Main.MesonDensity._strict_bw_stage_symbol(stage)
    density = if stage_sym === :strict_bw_stage1_reduced
        strict_bw_meson_density_summary(
            m_pi,
            gamma_pi,
            m_K,
            gamma_K,
            T_fm;
            μ_pi=Float64(μ_pi),
            μ_K=Float64(μ_K),
            d_pi=Int(d_pi),
            d_K=Int(d_K),
            qmax=qmax,
            q_nodes=q_nodes,
            omega_min=omega_min,
            omega_max=omega_max,
            omega_nodes=omega_nodes,
            gamma_zero_tol=gamma_zero_tol,
        )
    else
        thermo_params_norm = normalize_thermo_params(thermo_params)
        quark_params_raw = _require_result_field(meson_point, :quark_params)
        quark_params_norm = ensure_quark_params_has_A(
            normalize_quark_params(quark_params_raw),
            thermo_params_norm,
        )
        strict_bw_qpole_density_summary(
            m_pi,
            gamma_pi,
            m_K,
            gamma_K,
            quark_params_norm,
            thermo_params_norm;
            pi_channel=pi_channel,
            k_channel=k_channel,
            μ_pi=Float64(μ_pi),
            μ_K=Float64(μ_K),
            d_pi=Int(d_pi),
            d_K=Int(d_K),
            qmax=qmax,
            q_nodes=q_nodes,
            omega_max=omega_max,
            omega_nodes=omega_nodes,
            gamma_zero_tol=gamma_zero_tol,
            solver_iterations=solver_iterations,
            pole_residual_norm_max=pole_residual_norm_max,
            pole_require_converged=pole_require_converged,
        )
    end

    return merge(density, (
        T_fm=T_fm,
        xi=xi,
        stage=stage_sym,
        pi_channel=pi_channel,
        k_channel=k_channel,
        m_pi=m_pi,
        m_K=m_K,
        gamma_pi=gamma_pi,
        gamma_K=gamma_K,
        μ_pi=Float64(μ_pi),
        μ_K=Float64(μ_K),
        d_pi=Int(d_pi),
        d_K=Int(d_K),
    ))
end

function solve_gap_and_strict_bw_meson_density_point(
    T_fm::Real,
    mu_fm::Real;
    density_kwargs::NamedTuple=(;),
    kwargs...,
)
    meson_point = solve_gap_and_meson_point(T_fm, mu_fm; kwargs...)
    density = solve_strict_bw_meson_density_from_meson_point(meson_point; density_kwargs...)
    return merge(meson_point, (strict_bw_meson_density=density,))
end

"""
    solve_phase_shift_meson_density_from_meson_point(meson_point; kwargs...) -> NamedTuple

消费 meson workflow 返回值，执行当前最小 Phase-E3 口径下的相移双积分介子数密度后处理。

当前约束：
- 仅支持 `xi = 0`
- 仅支持 `π/K` 聚合通道
- 积分方案固定为 GL + 硬截断
"""
function solve_phase_shift_meson_density_from_meson_point(
    meson_point;
    pi_channel::Symbol=:pi,
    k_channel::Symbol=:K,
    μ_pi::Float64=0.0,
    μ_K::Float64=0.0,
    d_pi::Integer=meson_degeneracy(:pi),
    d_K::Integer=meson_degeneracy(:K),
    scheme::Symbol=:current,
    qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    omega_min::Float64=0.05,
    omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    eta::Float64=1e-6,
    real_axis_mode::Symbol=:finite_eta,
    phase_convention::Symbol=:arg_propagator,
    phase_display::Symbol=:unwrapped,
    density_policy::Symbol=:strict_normal_domain,
    bose_x_min::Float64=0.0,
    noanom_policy::Symbol=:none,
)
    thermo_params_raw = _require_result_field(meson_point, :thermo_params)
    quark_params_raw = _require_result_field(meson_point, :quark_params)

    thermo_params = normalize_thermo_params(thermo_params_raw)
    quark_params = ensure_quark_params_has_A(
        normalize_quark_params(quark_params_raw),
        thermo_params,
    )

    density = phase_shift_meson_density_summary(
        quark_params,
        thermo_params;
        pi_channel=pi_channel,
        k_channel=k_channel,
        μ_pi=μ_pi,
        μ_K=μ_K,
        d_pi=Int(d_pi),
        d_K=Int(d_K),
        scheme=scheme,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_min,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        eta=eta,
        real_axis_mode=real_axis_mode,
        phase_convention=phase_convention,
        phase_display=phase_display,
        density_policy=density_policy,
        bose_x_min=bose_x_min,
        noanom_policy=noanom_policy,
    )

    return merge(density, (
        m_pi=haskey(meson_point.meson_results, pi_channel) ? Float64(meson_point.meson_results[pi_channel].mass) : NaN,
        m_K=haskey(meson_point.meson_results, k_channel) ? Float64(meson_point.meson_results[k_channel].mass) : NaN,
    ))
end

function solve_gap_and_phase_shift_meson_density_point(
    T_fm::Real,
    mu_fm::Real;
    density_kwargs::NamedTuple=(;),
    kwargs...,
)
    meson_point = solve_gap_and_meson_point(T_fm, mu_fm; kwargs...)
    density = solve_phase_shift_meson_density_from_meson_point(meson_point; density_kwargs...)
    return merge(meson_point, (phase_shift_meson_density=density,))
end

function solve_phase_shift_derivative_reference_from_meson_point(
    meson_point;
    pi_channel::Symbol=:pi,
    k_channel::Symbol=:K,
    μ_pi::Float64=0.0,
    μ_K::Float64=0.0,
    d_pi::Integer=meson_degeneracy(:pi),
    d_K::Integer=meson_degeneracy(:K),
    scheme::Symbol=:current,
    qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    omega_min::Float64=0.05,
    omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    eta::Float64=1e-6,
    real_axis_mode::Symbol=:finite_eta,
)
    thermo_params_raw = _require_result_field(meson_point, :thermo_params)
    quark_params_raw = _require_result_field(meson_point, :quark_params)

    thermo_params = normalize_thermo_params(thermo_params_raw)
    quark_params = ensure_quark_params_has_A(
        normalize_quark_params(quark_params_raw),
        thermo_params,
    )

    density = phase_shift_meson_density_derivative_reference_summary(
        quark_params,
        thermo_params;
        pi_channel=pi_channel,
        k_channel=k_channel,
        μ_pi=μ_pi,
        μ_K=μ_K,
        d_pi=Int(d_pi),
        d_K=Int(d_K),
        scheme=scheme,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_min,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        eta=eta,
        real_axis_mode=real_axis_mode,
    )

    return merge(density, (
        m_pi=haskey(meson_point.meson_results, :pi) ? Float64(meson_point.meson_results[:pi].mass) : NaN,
        m_K=haskey(meson_point.meson_results, :K) ? Float64(meson_point.meson_results[:K].mass) : NaN,
    ))
end

function solve_phase_shift_point_diagnostic_from_meson_point(
    meson_point;
    mesons::Tuple=(:pi, :K),
    q_values::AbstractVector{<:Real}=[0.0],
    omega_values::AbstractVector{<:Real}=[0.2],
    scheme::Symbol=:current,
    eta::Float64=1e-6,
    real_axis_mode::Symbol=:finite_eta,
    phase_convention::Symbol=:arg_propagator,
    fd_step::Float64=1e-5,
)
    thermo_params_raw = _require_result_field(meson_point, :thermo_params)
    quark_params_raw = _require_result_field(meson_point, :quark_params)

    thermo_params = normalize_thermo_params(thermo_params_raw)
    quark_params = ensure_quark_params_has_A(
        normalize_quark_params(quark_params_raw),
        thermo_params,
    )

    rows = NamedTuple[]
    for meson in mesons, q in q_values, ω in omega_values
        push!(rows, phase_shift_point_diagnostic(
            meson,
            Float64(ω),
            Float64(q),
            quark_params,
            thermo_params;
            scheme=scheme,
            eta=eta,
            real_axis_mode=real_axis_mode,
            phase_convention=phase_convention,
            fd_step=fd_step,
        ))
    end

    return (
        T_fm=Float64(thermo_params.T),
        xi=Float64(thermo_params.ξ),
        scheme=scheme,
        eta=eta,
        real_axis_mode=real_axis_mode,
        phase_convention=phase_convention,
        fd_step=fd_step,
        rows=rows,
    )
end

end # module
