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

# ============================================================================
# PNJL 多种子初值（来自 legacy SeedStrategies.jl，用于 multi-seed 求解）
# ============================================================================

# 强子相典型种子（T<150 MeV 区间，禁闭相，Polyakov loop 接近零）
const _HADRON_SEED_5 = [-1.84329, -1.84329, -2.22701, 1.0e-5, 4.0e-5]
# 高温初值（T≈200 MeV，Polyakov loop ≈ 0.6）
const _HIGH_TEMP_SEED_5 = [-0.73192, -0.73192, -1.79539, 0.60532, 0.60532]
# 弱手征+禁闭（避免卡到坏分支）
const _WEAK_CHIRAL_CONF_SEED_5 = [-0.50, -0.50, -1.20, 1e-3, 1e-3]
# 人工高温候选
const _HT_GUESS_0p8_SEED_5 = [-0.50, -0.50, -1.20, 0.80, 0.80]
const _HT_GUESS_0p9_SEED_5 = [-0.30, -0.30, -0.90, 0.90, 0.90]
const _HT_GUESS_0p95_SEED_5 = [-0.20, -0.20, -0.70, 0.95, 0.95]

"""默认 PNJL 多种子列表（与 legacy MultiSeed() 一致）"""
const _PNJL_MULTI_SEEDS = [
    _HADRON_SEED_5,
    _HIGH_TEMP_SEED_5,
    _WEAK_CHIRAL_CONF_SEED_5,
    _HT_GUESS_0p8_SEED_5,
    _HT_GUESS_0p9_SEED_5,
    _HT_GUESS_0p95_SEED_5,
]

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

@inline gap_state_dim(::NJL2Model) = 2

@inline gap_state_dim(::AbstractPNJLModel) = 5

"""gap_initial_guess(model, T, mu_vec) -> AbstractVector

返回默认初值（用于 NLsolve）。

说明：这里给一个偏"强子相"的经验种子（与 legacy PNJL 种子尺度一致）。
"""
@inline function gap_initial_guess(::AbstractNJLModel, T, mu_vec)
    _ = (T, mu_vec)
    return [-1.84329, -1.84329, -2.22701]
end

@inline function gap_initial_guess(::NJL2Model, T, mu_vec)
    _ = (T, mu_vec)
    return [-1.84329, -1.84329]
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

@inline function _omega_worldsafe(model, st, T, μ; p_num::Int, t_num::Int, xi, kwargs...)
    return omega(model, st, T, μ; p_num=p_num, t_num=t_num, xi=xi, kwargs...)
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

    omega_fn = y -> begin
        Ty = eltype(y)
        if N == 3
            st = MeanFieldState(SVector{3, Ty}(y[1], y[2], y[3]); Phi=one(Ty), PhiBar=one(Ty))
            return _omega_worldsafe(model, st, T, μ; p_num=p_num, t_num=t_num, xi=xi, kwargs...)
        else
            st = MeanFieldState(SVector{3, Ty}(y[1], y[2], y[3]); Phi=y[4], PhiBar=y[5])
            return _omega_worldsafe(model, st, T, μ; p_num=p_num, t_num=t_num, xi=xi, kwargs...)
        end
    end

    g = ForwardDiff.gradient(omega_fn, x_state)
    return SVector{N, eltype(g)}(Tuple(g))
end

function gap_residual(
    model::NJL2Model,
    x,
    T,
    mu_vec;
    p_num::Int=64,
    t_num::Int=8,
    xi=0.0,
    kwargs...
)
    μ = normalize_mu_vec(mu_vec)

    x_state = if x isa MeanFieldState
        v = state_vector(x)
        SVector{2, eltype(v)}(v[1], v[2])
    else
        length(x) >= 2 || throw(ArgumentError("x must have length >= 2, got $(length(x))"))
        Tx = typeof(x[1])
        SVector{2, Tx}(x[1], x[2])
    end

    omega_fn = y -> begin
        Ty = eltype(y)
        st = MeanFieldState(SVector{3, Ty}(y[1], y[2], zero(Ty)); Phi=one(Ty), PhiBar=one(Ty))
        return omega(model, st, T, μ; p_num=p_num, t_num=t_num, xi=xi, kwargs...)
    end

    g = ForwardDiff.gradient(omega_fn, x_state)
    return SVector{2, eltype(g)}(Tuple(g))
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

@inline function _as_seed_vec2(seed)
    if seed isa MeanFieldState
        return [Float64(seed.phi[1]), Float64(seed.phi[2])]
    end
    if seed isa AbstractVector
        length(seed) >= 2 || throw(ArgumentError("initial_guess must have length >= 2"))
        return [Float64(seed[1]), Float64(seed[2])]
    end
    if seed isa NamedTuple
        st = MeanFieldState(seed)
        return [Float64(st.phi[1]), Float64(st.phi[2])]
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
    model::NJL2Model,
    T,
    mu_vec;
    solver::AbstractGapSolver=NLsolveGapSolver(),
    initial_guess=nothing,
    residual_norm_max::Real=1e-6,
    kwargs...
)
    x0 = initial_guess === nothing ? gap_initial_guess(model, T, mu_vec) : _as_seed_vec2(initial_guess)

    function residual!(F, x)
        r = gap_residual(model, x, T, mu_vec; kwargs...)
        @inbounds for i in 1:2
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
    phi = SVector{3, Float64}(x_sol[1], x_sol[2], 0.0)
    return MeanFieldState(phi; Phi=1.0, PhiBar=1.0)
end

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
- 当未提供 initial_guess 时，使用多种子策略（与 legacy MultiSeed 一致）：
  尝试多个初始猜测（含 Newton→Trust-Region 兜底），选择 Ω 最小的物理解。
- 该方法作为"可用的通用实现"，更具体模型（例如 `PNJLModel`）可继续提供更专用的求解器。
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

    s = solver isa NLsolveGapSolver ? solver : NLsolveGapSolver()

    function residual!(F, x)
        r = gap_residual(model, x, T, mu_vec; kwargs...)
        @inbounds for i in 1:5
            F[i] = convert(eltype(F), r[i])
        end
        return nothing
    end

    # 如果提供了 initial_guess，走单种子路径（与旧行为一致）
    if initial_guess !== nothing
        x0 = _as_seed_vec5(initial_guess)
        return _solve_gap_single_seed_pnjl(model, T, mu_vec, residual!, x0, s, residual_norm_max; kwargs...)
    end

    # 多种子路径：尝试所有默认种子，收集收敛的物理解，选 Ω 最小者
    best_state = nothing
    best_omega = Inf
    any_converged = false

    for seed in _PNJL_MULTI_SEEDS
        x0 = copy(seed)
        st = _try_solve_gap_pnjl(model, T, mu_vec, residual!, x0, s, residual_norm_max; kwargs...)
        st === nothing && continue

        any_converged = true
        ω = _safe_omega(model, st, T, mu_vec; kwargs...)
        if ω < best_omega
            best_omega = ω
            best_state = st
        end
    end

    if !any_converged
        error("solve_gap did not converge: all $(length(_PNJL_MULTI_SEEDS)) seeds failed at T=$T")
    end

    return best_state
end

"""单种子求解（无兜底，失败直接报错）"""
function _solve_gap_single_seed_pnjl(model, T, mu_vec, residual!, x0, s, residual_norm_max; kwargs...)
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

"""尝试用一个种子求解，Newton 失败则兜底 Trust-Region。返回 MeanFieldState 或 nothing。"""
function _try_solve_gap_pnjl(model, T, mu_vec, residual!, x0, s, residual_norm_max; kwargs...)
    autodiff_mode = s.jacobian === :forward ? :forward : :finite

    # Newton 尝试
    st = _try_nlsolve_pnjl(residual!, x0, autodiff_mode, s.method, s.xtol, s.ftol, residual_norm_max)
    if st !== nothing && _is_physical_pnjl(st)
        return st
    end

    # Trust-Region 兜底
    if s.method !== :trust_region
        st_tr = _try_nlsolve_pnjl(residual!, copy(x0), autodiff_mode, :trust_region, s.xtol, s.ftol, residual_norm_max)
        if st_tr !== nothing && _is_physical_pnjl(st_tr)
            return st_tr
        end
        # 如果 trust-region 收敛但不物理，而 Newton 收敛了（即使不物理），返回 Newton 结果
        if st !== nothing
            return st
        end
        if st_tr !== nothing
            return st_tr
        end
    end

    return st  # 可能是 nothing
end

"""单次 nlsolve 尝试，成功返回 MeanFieldState，失败返回 nothing"""
function _try_nlsolve_pnjl(residual!, x0, autodiff_mode, method, xtol, ftol, residual_norm_max)
    try
        result = nlsolve(residual!, x0; autodiff=autodiff_mode, method=method, xtol=xtol, ftol=ftol)
        resid = Float64(result.residual_norm)
        if result.f_converged && isfinite(resid) && resid <= Float64(residual_norm_max)
            x_sol = result.zero
            phi = SVector{3, Float64}(x_sol[1], x_sol[2], x_sol[3])
            return MeanFieldState(phi; Phi=Float64(x_sol[4]), PhiBar=Float64(x_sol[5]))
        end
    catch
    end
    return nothing
end

"""检查 PNJL 解的物理性：Φ, Φbar ∈ [0,1]，质量有限正定。"""
function _is_physical_pnjl(st::MeanFieldState)
    Φ  = st.Phi
    Φb = st.PhiBar
    tol = 0.01
    (Φ >= -tol && Φ <= 1.0 + tol) || return false
    (Φb >= -tol && Φb <= 1.0 + tol) || return false
    # 质量检查（粗略）：有效质量应为正 — 但此处只有 φ，不做质量计算
    all(isfinite, st.phi) || return false
    return true
end

"""安全计算 Ω，失败返回 Inf（便于比较）"""
function _safe_omega(model, st::MeanFieldState, T, mu_vec; kwargs...)
    try
        x_sv = state_vector(st)
        μ = normalize_mu_vec(mu_vec)
        ω = omega(model, x_sv, T, μ; kwargs...)
        return isfinite(ω) ? ω : Inf
    catch
        return Inf
    end
end
