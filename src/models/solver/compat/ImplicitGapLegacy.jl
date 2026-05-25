"""implicit_gap.jl

Models 侧的 legacy ImplicitDifferentiation 兼容层。

目标：在不对 NLsolve 迭代过程求导的前提下，通过隐函数定理获取 gap 解
x(T, μ) 的导数（例如 dφ/dT, dφ/dμ）。

设计对齐 legacy：
- legacy PNJL: Conditions.gap_conditions 使用 ForwardDiff.gradient 构造 F(x;θ)=0，
  NLsolve 使用 autodiff=:forward（等价于 Hessian 结构）。
- models: 通过 Conditions.gap_core_residual! 统一构建条件函数（NJL2/NJL3/PNJL 共用入口）。

注意：
- forward_solve_impl 只用于求“数值解”（primal），因此将 θ 强制转为 Float64 是 OK 的。
- conditions_impl 必须保持对 Dual 友好，避免 Float64(...) / float(...) 强制转换。
"""

using StaticArrays

export solve_pnjl_with_derivatives
export solve_pnjl_with_flavor_mu_derivatives
export derive_vec, derive_named

# Compat-only guard:
# This file is a legacy compatibility bridge and must not become a new primary
# solver/diff implementation surface.
const IMPLICIT_GAP_LEGACY_COMPAT_ONLY = true
const _IMPLICIT_COMPAT_DERIVATIVE_BACKENDS = (:auto, :taylordiff)

@inline function _validate_implicit_compat_derivative_backend(derivative_backend::Symbol)
    if derivative_backend === :forwarddiff
        throw(ArgumentError("derivative_backend=:forwarddiff has been retired from PNJL derivative wrappers; use derivative_backend=:auto or :taylordiff. Low-level create_*_implicit* factories remain qualified compat-only entrypoints for explicit legacy reference tests."))
    end
    derivative_backend in _IMPLICIT_COMPAT_DERIVATIVE_BACKENDS && return derivative_backend
    throw(ArgumentError("derivative_backend must be one of $(_IMPLICIT_COMPAT_DERIVATIVE_BACKENDS), got $(derivative_backend)"))
end

@inline _implicit_compat_backend(derivative_backend::Symbol) =
    _validate_implicit_compat_derivative_backend(derivative_backend) === :auto ? :taylordiff : derivative_backend

@inline function _reject_unsupported_td_wrapper_kwargs(kwargs)
    isempty(kwargs) && return nothing
    names = join(string.(keys(kwargs)), ", ")
    throw(ArgumentError("unsupported keyword(s) for TD-only PNJL derivative wrapper: $names. Legacy solver-specific keywords now apply only to qualified compat factories such as Models.create_implicit_gap_solver."))
end

@inline function _legacy_adapter_model_kind(model::AbstractQCDModel)
    return model isa RPNJLModel ? :RPNJL : :PNJL
end

@inline function symmetric_mu_direction_derivative(dx_dmu_vec::AbstractMatrix)
    size(dx_dmu_vec, 2) == 3 || throw(ArgumentError("dx_dmu_vec must have 3 columns for (μ_u, μ_d, μ_s), got size=$(size(dx_dmu_vec))"))
    return vec(sum(dx_dmu_vec; dims=2))
end

@inline function _normalize_flavor_mu_vec(mu_vec)
    μ = normalize_mu_vec(mu_vec)
    return [Float64(μ[1]), Float64(μ[2]), Float64(μ[3])]
end

@inline _td_zero_mu_direction() = SVector{3, Float64}(0.0, 0.0, 0.0)
@inline _td_symmetric_mu_direction() = SVector{3, Float64}(1.0, 1.0, 1.0)

@inline function _td_series_derivative_vector(v::SVector{N}, order::Int) where {N}
    return collect(SVector{N, Float64}(ntuple(i -> PNJLChiBTaylorDiff.nth_derivative_from_series(v[i], order), Val(N))))
end

function _td_gap_state_series(
    model::AbstractPNJLModel,
    T_fm::Real,
    mu_vec0,
    T_direction::Real,
    mu_direction;
    order::Int,
    xi::Real,
    p_num::Int,
    t_num::Int,
    series_iterations::Union{Nothing, Int},
    linear_solve::Symbol,
    series_residual_tol::Real,
)
    return PNJLChiBTaylorDiff.gap_series_parameter_direction(
        T_fm,
        mu_vec0,
        T_direction,
        mu_direction;
        order=order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=model,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
end

function _solve_pnjl_with_derivatives_taylordiff(
    model::AbstractPNJLModel,
    T_fm::Real,
    μ_fm::Real;
    order::Int,
    xi::Real,
    p_num::Int,
    t_num::Int,
    series_iterations::Union{Nothing, Int},
    linear_solve::Symbol,
    series_residual_tol::Real,
)
    order == 1 || order == 2 || throw(ArgumentError("order must be 1 or 2, got $order"))
    mu_vec0 = SVector{3, Float64}(Float64(μ_fm), Float64(μ_fm), Float64(μ_fm))

    s_T = _td_gap_state_series(
        model, T_fm, mu_vec0, 1.0, _td_zero_mu_direction();
        order=order, xi=xi, p_num=p_num, t_num=t_num,
        series_iterations=series_iterations, linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
    s_μ = _td_gap_state_series(
        model, T_fm, mu_vec0, 0.0, _td_symmetric_mu_direction();
        order=order, xi=xi, p_num=p_num, t_num=t_num,
        series_iterations=series_iterations, linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )

    x = _td_series_derivative_vector(s_T.x_state, 0)
    dx_dT = _td_series_derivative_vector(s_T.x_state, 1)
    dx_dμ = _td_series_derivative_vector(s_μ.x_state, 1)

    if order == 1
        return (x=x, dx_dT=dx_dT, dx_dμ=dx_dμ)
    end

    s_Tμ = _td_gap_state_series(
        model, T_fm, mu_vec0, 1.0, _td_symmetric_mu_direction();
        order=2, xi=xi, p_num=p_num, t_num=t_num,
        series_iterations=series_iterations, linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )

    d2x_dT2 = _td_series_derivative_vector(s_T.x_state, 2)
    d2x_dμ2 = _td_series_derivative_vector(s_μ.x_state, 2)
    d2x_Tμ_total = _td_series_derivative_vector(s_Tμ.x_state, 2)
    d2x_dTdμ = (d2x_Tμ_total .- d2x_dT2 .- d2x_dμ2) ./ 2

    return (x=x, dx_dT=dx_dT, dx_dμ=dx_dμ,
            d2x_dT2=d2x_dT2, d2x_dμ2=d2x_dμ2, d2x_dTdμ=d2x_dTdμ)
end

function _solve_pnjl_with_flavor_mu_derivatives_taylordiff(
    model::AbstractPNJLModel,
    T_fm::Real,
    mu_vec;
    order::Int,
    xi::Real,
    p_num::Int,
    t_num::Int,
    series_iterations::Union{Nothing, Int},
    linear_solve::Symbol,
    series_residual_tol::Real,
)
    order == 1 || throw(ArgumentError("solve_pnjl_with_flavor_mu_derivatives currently supports order=1 only"))
    μ0_arr = _normalize_flavor_mu_vec(mu_vec)
    μ0 = SVector{3, Float64}(μ0_arr[1], μ0_arr[2], μ0_arr[3])

    s_T = _td_gap_state_series(
        model, T_fm, μ0, 1.0, _td_zero_mu_direction();
        order=1, xi=xi, p_num=p_num, t_num=t_num,
        series_iterations=series_iterations, linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
    x = _td_series_derivative_vector(s_T.x_state, 0)
    dx_dT = _td_series_derivative_vector(s_T.x_state, 1)

    dx_dmu_vec = Matrix{Float64}(undef, 5, 3)
    for j in 1:3
        direction = SVector{3, Float64}(ntuple(i -> i == j ? 1.0 : 0.0, Val(3)))
        s_μj = _td_gap_state_series(
            model, T_fm, μ0, 0.0, direction;
            order=1, xi=xi, p_num=p_num, t_num=t_num,
            series_iterations=series_iterations, linear_solve=linear_solve,
            series_residual_tol=series_residual_tol,
        )
        dx_dmu_vec[:, j] = _td_series_derivative_vector(s_μj.x_state, 1)
    end

    return (x=x, mu_vec=μ0, dx_dT=dx_dT, dx_dmu_vec=dx_dmu_vec)
end

function build_pnjl_fixedmu_adapters(
    model::AbstractPNJLModel;
    xi::Real=0.0,
    p_num::Int=64,
    t_num::Int=8,
    kwargs...
)
    thermal_nodes = cached_nodes(p_num, t_num; p_max_inv_fm=Models.thermal_p_max_inv_fm(model))

    forward_solve = function (θ::AbstractVector)
        T_fm = Float64(θ[1])
        μ_fm = Float64(θ[2])
        mu_vec = SVector{3, Float64}(μ_fm, μ_fm, μ_fm)

        st = solve_gap(model, T_fm, mu_vec; xi=xi, p_num=p_num, t_num=t_num, kwargs...)
        return (collect(state_vector(st)), nothing)
    end

    conditions = function (θ::AbstractVector, x::AbstractVector, z)
        _ = z
        T_fm = θ[1]
        μ_fm = θ[2]
        mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)
        Tx = typeof(x[1])
        x_state = SVector{5, Tx}(Tuple(x))
        params = GapParams(T_fm, thermal_nodes, xi;
            p_num=p_num,
            t_num=t_num,
            model_kind=_legacy_adapter_model_kind(model),
        )

        Tout = promote_type(Tx, typeof(T_fm), typeof(mu_vec[1]))
        out = Vector{Tout}(undef, 5)
        gap_core_residual!(out, model, x_state, mu_vec, params)
        return out
    end

    return (forward_solve=forward_solve, conditions=conditions)
end

function build_pnjl_flavor_mu_adapters(
    model::AbstractPNJLModel;
    xi::Real=0.0,
    p_num::Int=64,
    t_num::Int=8,
    kwargs...
)
    thermal_nodes = cached_nodes(p_num, t_num; p_max_inv_fm=Models.thermal_p_max_inv_fm(model))

    forward_solve = function (θ::AbstractVector)
        length(θ) == 4 || throw(ArgumentError("flavor-mu implicit solver expects θ=[T, μ_u, μ_d, μ_s]"))
        T_fm = Float64(θ[1])
        mu_vec = SVector{3, Float64}(Float64(θ[2]), Float64(θ[3]), Float64(θ[4]))

        st = solve_gap(model, T_fm, mu_vec; xi=xi, p_num=p_num, t_num=t_num, kwargs...)
        return (collect(state_vector(st)), nothing)
    end

    conditions = function (θ::AbstractVector, x::AbstractVector, z)
        _ = z
        length(θ) == 4 || throw(ArgumentError("flavor-mu implicit solver expects θ=[T, μ_u, μ_d, μ_s]"))
        T_fm = θ[1]
        mu_vec = SVector{3}(θ[2], θ[3], θ[4])
        Tx = typeof(x[1])
        x_state = SVector{5, Tx}(Tuple(x))
        params = GapParams(T_fm, thermal_nodes, xi;
            p_num=p_num,
            t_num=t_num,
            model_kind=_legacy_adapter_model_kind(model),
        )

        Tout = promote_type(Tx, typeof(T_fm), typeof(mu_vec[1]))
        out = Vector{Tout}(undef, 5)
        gap_core_residual!(out, model, x_state, mu_vec, params)
        return out
    end

    return (forward_solve=forward_solve, conditions=conditions)
end

"""create_implicit_gap_solver(model; kwargs...) -> ImplicitFunction

Compat-only: 创建一个隐函数求解器：θ=[T, μ] -> x=[φu, φd, φs]。

该工厂已从公共 export 面降级，只供显式 legacy reference 测试或诊断脚本
qualified 调用。PNJL 生产导数入口请使用 `solve_pnjl_with_derivatives`
默认 TD 路线。

关键 kwargs：
- `xi`, `p_num`, `t_num`: 传给 omega/gap_residual
- `solver`: 传给 solve_gap（仅用于 primal forward solve）
- 其它 kwargs：也会透传给 solve_gap 与 gap_residual（如 residual_norm_max 等）
"""
function create_implicit_gap_solver(
    model::NJL2Model;
    xi::Real=0.0,
    p_num::Int=64,
    t_num::Int=8,
    solver::AbstractGapSolver=NLsolveGapSolver(),
    kwargs...
)
    gap_state_dim(model) == 2 || throw(ArgumentError("create_implicit_gap_solver(::NJL2Model) expects dim=2"))
    problem = build_njl_problem(
        model;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        solver=solver,
        kwargs...,
    )
    return build_implicit_solver(problem)
end

@inline function _pnjl_model_kind(thermo_backend::Symbol)
    if thermo_backend === :legacy
        return :LegacyPNJL
    elseif thermo_backend === :models
        return :PNJL
    end
    throw(ArgumentError("unknown thermo_backend=$thermo_backend (expected :legacy or :models)"))
end

"""create_pnjl_implicit_solver(; kwargs...) -> ImplicitFunction

Compat-only: 基于 models 入口创建 PNJL 5 维隐函数求解器：
- 参数 θ = [T, μ]
- 解向量 x = [φu, φd, φs, Φ, Φbar]

默认使用 models backend；可通过 `thermo_backend/solver_backend` 切换。
该入口不再作为推荐公开导数入口导出。
"""
function create_pnjl_implicit_solver(;
    xi::Real=0.0,
    p_num::Int=64,
    t_num::Int=8,
    thermo_backend::Symbol=:models,
    solver_backend::Symbol=:models,
    kwargs...
)
    kind = thermo_backend === :legacy ? :PNJL : _pnjl_model_kind(thermo_backend)
    model = create_model(kind)
    return create_implicit_gap_solver(
        model;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        solver_backend=solver_backend,
        kwargs...,
    )
end

"""solve_pnjl_with_derivatives(T_fm, μ_fm; kwargs...) -> NamedTuple

通过 models 求解器计算 PNJL 解及其导数。使用 TaylorDiff explicit
Taylor-series gap Newton；旧 `derivative_backend=:forwarddiff`
ImplicitFunction fallback 已下线：
- `order=1`：返回 `x, dx_dT, dx_dμ`
- `order=2`：额外返回 `d2x_dT2, d2x_dμ2, d2x_dTdμ`
"""
function solve_pnjl_with_derivatives(
    T_fm::Real,
    μ_fm::Real;
    order::Int=1,
    xi::Real=0.0,
    p_num::Int=64,
    t_num::Int=8,
    thermo_backend::Symbol=:models,
    solver_backend::Symbol=:models,
    derivative_backend::Symbol=:auto,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
    kwargs...
)
    backend = _implicit_compat_backend(derivative_backend)
    backend === :taylordiff || error("unreachable derivative backend: $backend")
    _reject_unsupported_td_wrapper_kwargs(kwargs)
    _ = solver_backend
    kind = thermo_backend === :legacy ? :PNJL : _pnjl_model_kind(thermo_backend)
    model = create_model(kind)
    return _solve_pnjl_with_derivatives_taylordiff(
        model,
        T_fm,
        μ_fm;
        order=order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
end

function derive_vec(
    model::AbstractPNJLModel,
    theta_vec::AbstractVector{<:Real};
    order::Int=1,
    xi::Real=0.0,
    p_num::Int=64,
    t_num::Int=8,
    thermo_backend::Symbol=:models,
    solver_backend::Symbol=:models,
    kwargs...
)
    length(theta_vec) == 2 || throw(ArgumentError("derive_vec expects theta_vec length 2 ([T_fm, μ_fm]), got $(length(theta_vec))"))
    _ = model
    return solve_pnjl_with_derivatives(
        theta_vec[1],
        theta_vec[2];
        order=order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        thermo_backend=thermo_backend,
        solver_backend=solver_backend,
        kwargs...,
    )
end

function derive_named(
    model::AbstractPNJLModel,
    theta_named::NamedTuple;
    kwargs...
)
    haskey(theta_named, :T_fm) || throw(ArgumentError("derive_named requires :T_fm"))
    haskey(theta_named, :μ_fm) || throw(ArgumentError("derive_named requires :μ_fm"))
    theta_vec = [theta_named[:T_fm], theta_named[:μ_fm]]
    return derive_vec(model, theta_vec; kwargs...)
end

"""create_flavor_mu_implicit_gap_solver(model::AbstractPNJLModel; kwargs...) -> ImplicitFunction

Compat-only: 创建一个 flavor 化学势版本的隐函数求解器：
- 参数 `θ = [T, μ_u, μ_d, μ_s]`
- 解向量 `x = [φu, φd, φs, Φ, Φbar]`

实现约定：
- `forward_solve_impl` 仅负责 primal solve，可安全使用 `Float64` 转换。
- `conditions_impl` 必须对 Dual 友好，不能把 `θ` 中的化学势分量压回 `Float64`。
该工厂已从公共 export 面降级，只供显式 legacy reference 测试或诊断脚本
qualified 调用。
"""
function create_flavor_mu_implicit_gap_solver(
    model::AbstractPNJLModel;
    xi::Real=0.0,
    p_num::Int=64,
    t_num::Int=8,
    kwargs...
)
    gap_state_dim(model) == 5 || throw(ArgumentError("create_flavor_mu_implicit_gap_solver(model::AbstractPNJLModel) expects dim=5"))
    problem = build_pnjl_flavor_mu_problem(
        model;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        kwargs...,
    )
    return build_implicit_solver(problem)
end

"""solve_pnjl_with_flavor_mu_derivatives(T_fm, mu_vec; order=1, kwargs...) -> NamedTuple

通过 flavor 化学势版本计算 PNJL 解及其导数，默认使用 TaylorDiff
方向 series：
- `x`: 5 维平衡态向量
- `dx_dT`: 对温度的一阶导数
- `dx_dmu_vec`: 对 `(μ_u, μ_d, μ_s)` 的 5×3 Jacobian

当前仅支持 `order=1`。旧 `derivative_backend=:forwarddiff` ImplicitFunction
fallback 已下线；更高阶张量留待 susceptibility 层真正需要时再接入，
避免在底层接口过早固化不稳定的张量布局。
"""
function solve_pnjl_with_flavor_mu_derivatives(
    T_fm::Real,
    mu_vec;
    order::Int=1,
    xi::Real=0.0,
    p_num::Int=64,
    t_num::Int=8,
    thermo_backend::Symbol=:models,
    solver_backend::Symbol=:models,
    derivative_backend::Symbol=:auto,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
    kwargs...
)
    order == 1 || throw(ArgumentError("solve_pnjl_with_flavor_mu_derivatives currently supports order=1 only"))

    kind = thermo_backend === :legacy ? :PNJL : _pnjl_model_kind(thermo_backend)
    model = create_model(kind)

    backend = _implicit_compat_backend(derivative_backend)
    backend === :taylordiff || error("unreachable derivative backend: $backend")
    _reject_unsupported_td_wrapper_kwargs(kwargs)
    _ = solver_backend
    return _solve_pnjl_with_flavor_mu_derivatives_taylordiff(
        model,
        T_fm,
        mu_vec;
        order=order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
end

function create_implicit_gap_solver(
    model::AbstractNJLModel;
    xi::Real=0.0,
    p_num::Int=64,
    t_num::Int=8,
    solver::AbstractGapSolver=NLsolveGapSolver(),
    kwargs...
)
    gap_state_dim(model) == 3 || throw(ArgumentError("create_implicit_gap_solver currently supports dim=3 only"))
    problem = build_njl_problem(
        model;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        solver=solver,
        kwargs...,
    )
    return build_implicit_solver(problem)
end

"""create_implicit_gap_solver(model::AbstractPNJLModel; kwargs...) -> ImplicitFunction

Compat-only: 创建一个隐函数求解器：θ=[T, μ] -> x=[φu, φd, φs, Φ, Φbar]。

实现说明：
- forward_solve_impl：调用 `solve_gap(model, T, μ)` 得到 `MeanFieldState`，再展开为 5 维向量。
- conditions_impl：调用 `gap_residual(model, x, T, μ)`。

注意：
- 目前默认假设对称化学势（μu=μd=μs），与 legacy FixedMu / 当前 PNJLModel.solve_gap 的限制一致。
- 该工厂已从公共 export 面降级；PNJL 导数生产入口请使用
  `solve_pnjl_with_derivatives` 的 TD 默认路径。
"""
function create_implicit_gap_solver(
    model::AbstractPNJLModel;
    xi::Real=0.0,
    p_num::Int=64,
    t_num::Int=8,
    kwargs...
)
    gap_state_dim(model) == 5 || throw(ArgumentError("create_implicit_gap_solver(model::AbstractPNJLModel) expects dim=5"))
    problem = build_pnjl_fixedmu_problem(
        model;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        kwargs...,
    )
    return build_implicit_solver(problem)
end
