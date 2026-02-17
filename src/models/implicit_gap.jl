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

    forward_solve_impl = function (θ::AbstractVector)
        T_fm = Float64(θ[1])
        μ_fm = Float64(θ[2])
        mu_vec = SVector{3, Float64}(μ_fm, μ_fm, μ_fm)

        st = solve_gap(model, T_fm, mu_vec; solver=solver, xi=xi, p_num=p_num, t_num=t_num, kwargs...)
        return ([st.phi[1], st.phi[2]], nothing)
    end

    conditions_impl = function (θ::AbstractVector, x::AbstractVector, z)
        _ = z
        T_fm = θ[1]
        μ_fm = θ[2]
        mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)

        r = gap_residual(model, x, T_fm, mu_vec;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            kwargs...)
        return Vector(r)
    end

    return ImplicitFunction(
        forward_solve_impl,
        conditions_impl;
        linear_solver=DirectLinearSolver(),
        representation=MatrixRepresentation(),
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

    # primal forward solve (θ -> x)
    forward_solve_impl = function (θ::AbstractVector)
        T_fm = Float64(θ[1])
        μ_fm = Float64(θ[2])
        mu_vec = SVector{3, Float64}(μ_fm, μ_fm, μ_fm)

        st = solve_gap(model, T_fm, mu_vec; solver=solver, xi=xi, p_num=p_num, t_num=t_num, kwargs...)
        return (collect(st.phi), nothing)
    end

    # conditions for implicit differentiation (θ, x) -> F(x;θ)
    conditions_impl = function (θ::AbstractVector, x::AbstractVector, z)
        _ = z
        T_fm = θ[1]
        μ_fm = θ[2]
        mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)

        r = gap_residual(model, x, T_fm, mu_vec;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            kwargs...)
        return Vector(r)
    end

    return ImplicitFunction(
        forward_solve_impl,
        conditions_impl;
        linear_solver=DirectLinearSolver(),
        representation=MatrixRepresentation(),
    )
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

    forward_solve_impl = function (θ::AbstractVector)
        T_fm = Float64(θ[1])
        μ_fm = Float64(θ[2])
        mu_vec = SVector{3, Float64}(μ_fm, μ_fm, μ_fm)

        st = solve_gap(model, T_fm, mu_vec; xi=xi, p_num=p_num, t_num=t_num, kwargs...)
        return (collect(state_vector(st)), nothing)
    end

    conditions_impl = function (θ::AbstractVector, x::AbstractVector, z)
        _ = z
        T_fm = θ[1]
        μ_fm = θ[2]
        mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)

        r = gap_residual(model, x, T_fm, mu_vec;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            kwargs...)
        return Vector(r)
    end

    return ImplicitFunction(
        forward_solve_impl,
        conditions_impl;
        linear_solver=DirectLinearSolver(),
        representation=MatrixRepresentation(),
    )
end
