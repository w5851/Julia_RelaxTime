"""implicit_gap.jl

Models 侧的 ImplicitDifferentiation 接口（阶段 1 组件）。

目标：在不对 NLsolve 迭代过程求导的前提下，通过隐函数定理获取 gap 解
x(T, μ) 的导数（例如 dφ/dT, dφ/dμ）。

设计对齐 legacy：
- legacy PNJL: Conditions.gap_conditions 使用 ForwardDiff.gradient 构造 F(x;θ)=0，
  NLsolve 使用 autodiff=:forward（等价于 Hessian 结构）。
- models: 复用 gap_residual(model, x, T, mu_vec) 作为条件函数。

注意：
- forward_solve_impl 只用于求“数值解”（primal），因此将 θ 强制转为 Float64 是 OK 的。
- conditions_impl 必须保持对 Dual 友好，避免 Float64(...) / float(...) 强制转换。
"""

using StaticArrays
using ForwardDiff
using ImplicitDifferentiation

export create_implicit_gap_solver
export create_flavor_mu_implicit_gap_solver
export create_pnjl_implicit_solver
export solve_pnjl_with_derivatives
export solve_pnjl_with_flavor_mu_derivatives
export build_pnjl_fixedmu_adapters
export build_pnjl_flavor_mu_adapters
export derive_vec, derive_named

@inline function symmetric_mu_direction_derivative(dx_dmu_vec::AbstractMatrix)
    size(dx_dmu_vec, 2) == 3 || throw(ArgumentError("dx_dmu_vec must have 3 columns for (μ_u, μ_d, μ_s), got size=$(size(dx_dmu_vec))"))
    return vec(sum(dx_dmu_vec; dims=2))
end

@inline function _normalize_flavor_mu_vec(mu_vec)
    μ = normalize_mu_vec(mu_vec)
    return [Float64(μ[1]), Float64(μ[2]), Float64(μ[3])]
end

function build_pnjl_fixedmu_adapters(
    model::AbstractPNJLModel;
    xi::Real=0.0,
    p_num::Int=64,
    t_num::Int=8,
    kwargs...
)
    thermal_nodes = cached_nodes(p_num, t_num)

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
        params = GapParams(T_fm, thermal_nodes, xi; p_num=p_num, t_num=t_num, model_kind=:PNJL)

        Tout = promote_type(Tx, typeof(T_fm), typeof(mu_vec[1]))
        out = Vector{Tout}(undef, 5)
        gap_core_residual!(out, x_state, mu_vec, params)
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
    thermal_nodes = cached_nodes(p_num, t_num)

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
        params = GapParams(T_fm, thermal_nodes, xi; p_num=p_num, t_num=t_num, model_kind=:PNJL)

        Tout = promote_type(Tx, typeof(T_fm), typeof(mu_vec[1]))
        out = Vector{Tout}(undef, 5)
        gap_core_residual!(out, x_state, mu_vec, params)
        return out
    end

    return (forward_solve=forward_solve, conditions=conditions)
end

"""create_implicit_gap_solver(model; kwargs...) -> ImplicitFunction

创建一个隐函数求解器：θ=[T, μ] -> x=[φu, φd, φs]，并可通过 AD 获取 dx/dθ。

当前仅对 `AbstractNJLModel` 提供 3 维实现；PNJL 家族可在后续扩展到 5 维或更高。

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

基于 models 入口创建 PNJL 5 维隐函数求解器：
- 参数 θ = [T, μ]
- 解向量 x = [φu, φd, φs, Φ, Φbar]

默认使用 models backend；可通过 `thermo_backend/solver_backend` 切换。
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

通过 models 隐函数求解器计算 PNJL 解及其导数：
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
    kwargs...
)
    solver = create_pnjl_implicit_solver(
        ;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        thermo_backend=thermo_backend,
        solver_backend=solver_backend,
        kwargs...,
    )

    θ = [Float64(T_fm), Float64(μ_fm)]
    x, _ = solver(θ)

    if order == 1
        dx_dT = ForwardDiff.derivative(T -> solver([T, θ[2]])[1], θ[1])
        dx_dμ = ForwardDiff.derivative(μ -> solver([θ[1], μ])[1], θ[2])
        return (x=x, dx_dT=dx_dT, dx_dμ=dx_dμ)
    elseif order == 2
        dx_dT = ForwardDiff.derivative(T -> solver([T, θ[2]])[1], θ[1])
        dx_dμ = ForwardDiff.derivative(μ -> solver([θ[1], μ])[1], θ[2])

        d2x_dT2 = ForwardDiff.derivative(
            T -> ForwardDiff.derivative(t -> solver([t, θ[2]])[1], T),
            θ[1],
        )
        d2x_dμ2 = ForwardDiff.derivative(
            μ -> ForwardDiff.derivative(m -> solver([θ[1], m])[1], μ),
            θ[2],
        )
        d2x_dTdμ = ForwardDiff.derivative(
            T -> ForwardDiff.derivative(μ -> solver([T, μ])[1], θ[2]),
            θ[1],
        )

        return (x=x, dx_dT=dx_dT, dx_dμ=dx_dμ,
                d2x_dT2=d2x_dT2, d2x_dμ2=d2x_dμ2, d2x_dTdμ=d2x_dTdμ)
    else
        throw(ArgumentError("order must be 1 or 2, got $order"))
    end
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

创建一个 flavor 化学势版本的隐函数求解器：
- 参数 `θ = [T, μ_u, μ_d, μ_s]`
- 解向量 `x = [φu, φd, φs, Φ, Φbar]`

实现约定：
- `forward_solve_impl` 仅负责 primal solve，可安全使用 `Float64` 转换。
- `conditions_impl` 必须对 Dual 友好，不能把 `θ` 中的化学势分量压回 `Float64`。
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

通过 flavor 化学势版本的隐函数求解器计算 PNJL 解及其导数：
- `x`: 5 维平衡态向量
- `dx_dT`: 对温度的一阶导数
- `dx_dmu_vec`: 对 `(μ_u, μ_d, μ_s)` 的 5×3 Jacobian

当前仅支持 `order=1`。更高阶张量留待 susceptibility 层真正需要时再接入，
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
    kwargs...
)
    order == 1 || throw(ArgumentError("solve_pnjl_with_flavor_mu_derivatives currently supports order=1 only"))

    kind = thermo_backend === :legacy ? :PNJL : _pnjl_model_kind(thermo_backend)
    model = create_model(kind)

    solver = create_flavor_mu_implicit_gap_solver(
        model;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        solver_backend=solver_backend,
        kwargs...,
    )

    μ0 = _normalize_flavor_mu_vec(mu_vec)
    θ = [Float64(T_fm), μ0[1], μ0[2], μ0[3]]
    x, _ = solver(θ)

    dx_dT = ForwardDiff.derivative(T -> solver([T, μ0[1], μ0[2], μ0[3]])[1], θ[1])
    dx_dmu_vec = ForwardDiff.jacobian(μ -> solver([θ[1], μ[1], μ[2], μ[3]])[1], μ0)

    return (x=x, mu_vec=SVector{3}(μ0[1], μ0[2], μ0[3]), dx_dT=dx_dT, dx_dmu_vec=dx_dmu_vec)
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

创建一个隐函数求解器：θ=[T, μ] -> x=[φu, φd, φs, Φ, Φbar]。

实现说明：
- forward_solve_impl：调用 `solve_gap(model, T, μ)` 得到 `MeanFieldState`，再展开为 5 维向量。
- conditions_impl：调用 `gap_residual(model, x, T, μ)`。

注意：
- 目前默认假设对称化学势（μu=μd=μs），与 legacy FixedMu / 当前 PNJLModel.solve_gap 的限制一致。
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
