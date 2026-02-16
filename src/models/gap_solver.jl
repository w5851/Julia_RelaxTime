"""gap_solver.jl

阶段 0/1：models 子系统的通用 gap 求解入口。

设计目标：
- 求解器逻辑可复用（NLsolve/容差/初值/判据）
- 模型差异通过多重派发提供（残差/参数化/物理性检查）

当前实现策略（保守、可工作）：
- 对 NJL/PNJL：默认用 ForwardDiff 计算 ∇Ω=0（并保留有限差分兜底）
- 对 PNJL：仍保留在各自 `solve_gap(::PNJLModel, ...)` 中走 legacy 适配器
"""

using LinearAlgebra: norm
using NLsolve
using StaticArrays
using ForwardDiff

export NLsolveGapSolver

abstract type AbstractGapSolver end

"""NLsolveGapSolver

- `method`: `:newton` / `:trust_region` 等（传给 NLsolve）
- `xtol` / `ftol`: 传给 NLsolve
"""
Base.@kwdef struct NLsolveGapSolver <: AbstractGapSolver
    method::Symbol = :newton
    xtol::Float64 = 1e-10
    ftol::Float64 = 1e-10
    # Jacobian strategy for NLsolve.
    #
    # - `:forward`: let NLsolve compute ∂F/∂x via ForwardDiff.
    #   If F(x) is itself defined using ForwardDiff (e.g. F = ∇Ω), this implies
    #   nested Dual numbers and effectively gives a Hessian-like object.
    #   This is the same pattern used in legacy PNJL (gap_conditions uses a gradient,
    #   and NLsolve autodiff builds its Jacobian).
    #
    # - `:finite`: use finite-difference Jacobians inside NLsolve.
    #   This avoids nested AD and can be more robust, but is often slower.
    jacobian::Symbol = :forward
end

# -----------------------------------------------------------------------------
# Gap residual interface (per-model)
# -----------------------------------------------------------------------------

"""gap_state_dim(model) -> Int

返回 gap 求解的未知量维度。

阶段 0：
- NJL 默认 3（φ_u, φ_d, φ_s）
- PNJL/RPNJL 未来可扩展到 5（再加 Φ, Φbar）
"""
@inline gap_state_dim(::AbstractNJLModel) = 3

@inline gap_state_dim(::AbstractPNJLModel) = 5

"""gap_initial_guess(model, T, mu_vec) -> AbstractVector

返回默认初值（用于 NLsolve）。

说明：这里给一个偏“强子相”的经验种子（与 legacy PNJL 种子尺度一致）。
"""
@inline function gap_initial_guess(::AbstractNJLModel, T, mu_vec)
    _ = (T, mu_vec)
    return [-1.84329, -1.84329, -2.22701]
end

@inline function gap_initial_guess(::AbstractPNJLModel, T, mu_vec)
    _ = (T, mu_vec)
    # A conservative seed: NJL-like chiral condensates + mid-range Polyakov variables.
    return [-1.84329, -1.84329, -2.22701, 0.5, 0.5]
end

"""gap_residual(model, x, T, mu_vec; kwargs...) -> SVector

返回 gap 方程残差向量。

默认实现：通用驻点条件 ∇ₓΩ(x)=0。

- NJL（3 维）：x=(φ_u, φ_d, φ_s)，内部取 Φ=Φbar=1。
- PNJL（5 维）：x=(φ_u, φ_d, φ_s, Φ, Φbar)。

说明：
默认使用 ForwardDiff 计算梯度。

更具体模型可通过多重派发重载该方法。
"""
@inline function _as_svec(x, ::Val{N}) where {N}
    if x isa SVector{N}
        return x
    end
    Tx = typeof(x[1])
    if N == 3
        return SVector{3, Tx}(x[1], x[2], x[3])
    elseif N == 5
        return SVector{5, Tx}(x[1], x[2], x[3], x[4], x[5])
    end
    throw(ArgumentError("unsupported state dim N=$N"))
end

function gap_residual(
    model::AbstractQCDModel,
    x,
    T,
    mu_vec;
    p_num::Int=64,
    t_num::Int=8,
    xi=0.0,
    kwargs...
)
    μ = normalize_mu_vec(mu_vec)
    N = gap_state_dim(model)
    (N == 3 || N == 5) || throw(ArgumentError("gap_state_dim(model) must be 3 or 5, got $N"))

    x_state = if x isa MeanFieldState
        v = state_vector(x)
        if N == 3
            SVector{3, eltype(v)}(v[1], v[2], v[3])
        else
            v
        end
    else
        length(x) >= N || throw(ArgumentError("x must have length >= $N, got $(length(x))"))
        _as_svec(x, Val(N))
    end

    @inline function _omega_worldsafe(st)
        # Legacy adapters are loaded lazily; use invokelatest to avoid world-age
        # dispatch falling back to generic omega/omega_components during nested AD.
        if (isdefined(@__MODULE__, :LegacyPNJLModel) && model isa LegacyPNJLModel) ||
           (isdefined(@__MODULE__, :LegacyNJLModel) && model isa LegacyNJLModel)
            return Base.invokelatest(omega, model, st, T, μ; p_num=p_num, t_num=t_num, xi=xi, kwargs...)
        end
        return omega(model, st, T, μ; p_num=p_num, t_num=t_num, xi=xi, kwargs...)
    end

    omega_fn = y -> begin
        Ty = eltype(y)
        if N == 3
            st = MeanFieldState(SVector{3, Ty}(y[1], y[2], y[3]); Phi=one(Ty), PhiBar=one(Ty))
            return _omega_worldsafe(st)
        else
            st = MeanFieldState(SVector{3, Ty}(y[1], y[2], y[3]); Phi=y[4], PhiBar=y[5])
            return _omega_worldsafe(st)
        end
    end

    g = ForwardDiff.gradient(omega_fn, x_state)
    return SVector{N, eltype(g)}(Tuple(g))
end

# -----------------------------------------------------------------------------
# Generic solve_gap for models that implement gap_residual
# -----------------------------------------------------------------------------

@inline function _as_seed_vec3(seed)
    if seed isa MeanFieldState
        return Vector{Float64}(seed.phi)
    end
    if seed isa AbstractVector
        length(seed) >= 3 || throw(ArgumentError("initial_guess must have length >= 3"))
        return [Float64(seed[1]), Float64(seed[2]), Float64(seed[3])]
    end
    if seed isa NamedTuple
        st = MeanFieldState(seed)
        return Vector{Float64}(st.phi)
    end
    throw(ArgumentError("unsupported initial_guess type: $(typeof(seed))"))
end

@inline function _as_seed_vec5(seed)
    if seed isa MeanFieldState
        v = state_vector(seed)
        return [Float64(v[1]), Float64(v[2]), Float64(v[3]), Float64(v[4]), Float64(v[5])]
    end
    if seed isa AbstractVector
        length(seed) >= 5 || throw(ArgumentError("initial_guess must have length >= 5"))
        return [Float64(seed[1]), Float64(seed[2]), Float64(seed[3]), Float64(seed[4]), Float64(seed[5])]
    end
    if seed isa NamedTuple
        st = MeanFieldState(seed)
        v = state_vector(st)
        return [Float64(v[1]), Float64(v[2]), Float64(v[3]), Float64(v[4]), Float64(v[5])]
    end
    throw(ArgumentError("unsupported initial_guess type: $(typeof(seed))"))
end

"""solve_gap(model::AbstractNJLModel, T, mu_vec; kwargs...) -> MeanFieldState

通用求解入口：调用 `gap_residual(model, x, T, mu_vec; ...)` 并用 NLsolve 求根。

注意：PNJL 目前有更具体的 `solve_gap(::PNJLModel, ...)`（legacy 适配器），会优先走那个。
"""
function solve_gap(
    model::AbstractNJLModel,
    T,
    mu_vec;
    solver::AbstractGapSolver=NLsolveGapSolver(),
    initial_guess=nothing,
    residual_norm_max::Real=1e-6,
    kwargs...
)
    dim = gap_state_dim(model)
    dim == 3 || throw(ArgumentError("generic solve_gap currently supports dim=3 only, got $dim"))

    x0 = initial_guess === nothing ? gap_initial_guess(model, T, mu_vec) : _as_seed_vec3(initial_guess)

    # NLsolve expects an in-place residual! of Float64 arrays.
    function residual!(F, x)
        r = gap_residual(model, x, T, mu_vec; kwargs...)
        @inbounds for i in 1:3
            F[i] = convert(eltype(F), r[i])
        end
        return nothing
    end

    s = solver isa NLsolveGapSolver ? solver : NLsolveGapSolver()
    autodiff_mode = s.jacobian === :forward ? :forward : :finite
    result = nlsolve(residual!, x0; autodiff=autodiff_mode, method=s.method, xtol=s.xtol, ftol=s.ftol)

    resid = Float64(result.residual_norm)
    if !(result.f_converged && isfinite(resid) && resid <= Float64(residual_norm_max))
        error("solve_gap did not converge (residual_norm=$resid, f_converged=$(result.f_converged))")
    end

    x_sol = result.zero
    phi = SVector{3, Float64}(x_sol[1], x_sol[2], x_sol[3])
    return MeanFieldState(phi; Phi=1.0, PhiBar=1.0)
end

"""solve_gap(model::AbstractPNJLModel, T, mu_vec; kwargs...) -> MeanFieldState

通用 5 维 PNJL 家族求解入口：调用 `gap_residual(model, x, T, mu_vec; ...)` 并用 NLsolve 求根。

说明：
- 默认 residual 为 ∇Ω=0（AD 梯度），NLsolve 的 Jacobian 策略由 `solver.jacobian` 控制。
- 该方法作为“可用的通用实现”，更具体模型（例如 `PNJLModel`）可继续提供更专用的求解器。
"""
function solve_gap(
    model::AbstractPNJLModel,
    T,
    mu_vec;
    solver::AbstractGapSolver=NLsolveGapSolver(),
    initial_guess=nothing,
    residual_norm_max::Real=1e-6,
    kwargs...
)
    dim = gap_state_dim(model)
    dim == 5 || throw(ArgumentError("generic PNJL solve_gap expects dim=5, got $dim"))

    x0 = initial_guess === nothing ? gap_initial_guess(model, T, mu_vec) : _as_seed_vec5(initial_guess)

    function residual!(F, x)
        r = gap_residual(model, x, T, mu_vec; kwargs...)
        @inbounds for i in 1:5
            F[i] = convert(eltype(F), r[i])
        end
        return nothing
    end

    s = solver isa NLsolveGapSolver ? solver : NLsolveGapSolver()
    autodiff_mode = s.jacobian === :forward ? :forward : :finite
    result = nlsolve(residual!, x0; autodiff=autodiff_mode, method=s.method, xtol=s.xtol, ftol=s.ftol)

    resid = Float64(result.residual_norm)
    if !(result.f_converged && isfinite(resid) && resid <= Float64(residual_norm_max))
        error("solve_gap did not converge (residual_norm=$resid, f_converged=$(result.f_converged))")
    end

    x_sol = result.zero
    phi = SVector{3, Float64}(x_sol[1], x_sol[2], x_sol[3])
    return MeanFieldState(phi; Phi=Float64(x_sol[4]), PhiBar=Float64(x_sol[5]))
end
