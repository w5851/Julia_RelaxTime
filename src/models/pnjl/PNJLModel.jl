"""PNJLModel

新架构下的各向同性 PNJL 模型适配器。

当前阶段目标：开始“替换”而不大改旧实现 —— 直接复用旧路径
[src/pnjl/core/Thermodynamics.jl](../../pnjl/core/Thermodynamics.jl) 的公式实现，
对外提供新架构接口（供 `Models.omega` 组装）。

注意：这是过渡实现，后续可再把 Integrals/Polyakov 等下沉到 models/base。
"""

using StaticArrays
using Base.MathConstants: π

# Include-once helper
const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

# constants
const _CONSTANTS_PATH = normpath(joinpath(@__DIR__, "..", "..", "Constants_PNJL.jl"))
IncludeOnce.include_once!(Main, :Constants_PNJL, _CONSTANTS_PATH)
using Main.Constants_PNJL: pnjl_constants

include(joinpath(@__DIR__, "PNJLIntegrals.jl"))
const _pnjl_model_cached_nodes = PNJLIntegrals.cached_nodes
const _pnjl_model_log_sum = PNJLIntegrals.calculate_log_sum
const _PNJL_POLYAKOV_EPS = 1e-16

const _PNJL_SOLVER_PATH = normpath(joinpath(@__DIR__, "..", "..", "pnjl", "solver", "Solver.jl"))

@inline function _ensure_pnjl_solver_loaded()
    if !(isdefined(Main, :Solver) && isdefined(Main.Solver, :solve) && isdefined(Main.Solver, :FixedMu))
        Base.include(Main, _PNJL_SOLVER_PATH)
    end
    return nothing
end

export PNJLModel

"""各向同性色散 PNJL 模型（过渡期：参数来自 `Constants_PNJL.pnjl_constants`）。

阶段 2（配置注入）最小切片：
- 先把 models 侧的关键常量（如 `N_color`）从全局常量迁移到实例字段，避免依赖模块 load-time 缓存。
"""
struct PNJLModel <: AbstractPNJLModel
    consts::NamedTuple
end

function PNJLModel(; profile::String=get(ENV, "PNJL_PARAM_PROFILE", "default"), physics_profile::String=get(ENV, "PHYSICS_PARAM_PROFILE", "default"))
    return PNJLModel(pnjl_constants(profile=profile, physics_profile=physics_profile))
end

PNJLModel() = PNJLModel(; profile=get(ENV, "PNJL_PARAM_PROFILE", "default"), physics_profile=get(ENV, "PHYSICS_PARAM_PROFILE", "default"))

# -------------------------
# Interface implementations
# -------------------------

function calculate_mass_vec(model::PNJLModel, φ::SVector{3, T}; kwargs...) where {T}
    φ_u, φ_d, φ_s = φ
    m_ud0 = convert(T, Float64(model.consts.m_ud0_inv_fm))
    m_s0 = convert(T, Float64(model.consts.m_s0_inv_fm))
    G = convert(T, Float64(model.consts.G_fm2))
    K = convert(T, Float64(model.consts.K_fm5))

    return SVector{3, T}(
        m_ud0 - 4 * G * φ_u + 2 * K * φ_d * φ_s,
        m_ud0 - 4 * G * φ_d + 2 * K * φ_u * φ_s,
        m_s0 - 4 * G * φ_s + 2 * K * φ_u * φ_d,
    )
end

function calculate_chiral(model::PNJLModel, φ::SVector{3, T}; kwargs...) where {T}
    G = convert(T, Float64(model.consts.G_fm2))
    K = convert(T, Float64(model.consts.K_fm5))
    return 2 * G * sum(φ .^ 2) - 4 * K * prod(φ)
end

@inline function _pnjl_safe_log(x; min_val=_PNJL_POLYAKOV_EPS)
    min_x = one(x) * min_val
    x <= zero(x) && return log(min_x)
    return x < min_x ? log(min_x) : log(x)
end

function polyakov_potential(model::PNJLModel, Φ, Φbar, T_fm; kwargs...)
    TT = promote_type(typeof(Φ), typeof(Φbar), typeof(T_fm))
    ΦT = convert(TT, Φ)
    ΦbarT = convert(TT, Φbar)
    TT_fm = convert(TT, T_fm)

    T0 = convert(TT, Float64(model.consts.T0_inv_fm))
    a0 = convert(TT, Float64(model.consts.a0))
    a1 = convert(TT, Float64(model.consts.a1))
    a2 = convert(TT, Float64(model.consts.a2))
    b3 = convert(TT, Float64(model.consts.b3))

    T_ratio = T0 / TT_fm
    Ta = a0 + a1 * T_ratio + a2 * T_ratio^2
    Tb = b3 * T_ratio^3
    value = 1 - 6 * ΦbarT * ΦT + 4 * (ΦbarT^3 + ΦT^3) - 3 * (ΦbarT * ΦT)^2
    return TT_fm^4 * (-0.5 * Ta * ΦbarT * ΦT + Tb * _pnjl_safe_log(value))
end

@inline function _vacuum_integral_with_cutoff(mass::TF, Λ::TF) where {TF}
    mass_abs = abs(mass)
    epsilon = one(mass_abs) * 1e-12
    mass_safe = mass_abs + epsilon
    sqrt_term = sqrt(Λ^2 + mass_safe^2)
    poly_part = Λ * sqrt_term * (2 * Λ^2 + mass_safe^2)
    log_term = mass_safe^4 * log((Λ + sqrt_term) / mass_safe)
    return (poly_part - log_term) / (16 * π^2)
end

function vacuum_contribution(model::PNJLModel, masses::SVector{3, T}; kwargs...) where {T}
    Λ = convert(T, Float64(model.consts.Λ_inv_fm))
    Nc = convert(T, Int(model.consts.N_color))

    total = zero(T)
    @inbounds for i in 1:3
        total += _vacuum_integral_with_cutoff(masses[i], Λ)
    end
    return -2 * Nc * total
end

function thermal_contribution(
    ::PNJLModel,
    masses::SVector{3, T},
    Φ,
    Φbar,
    mu_vec,
    T_fm;
    p_num::Int=PNJLIntegrals.DEFAULT_MOMENTUM_COUNT,
    t_num::Int=PNJLIntegrals.DEFAULT_THETA_COUNT,
    xi=0.0,
    kwargs...
) where {T}
    p_mesh, cosθ_mesh, coefficients = _pnjl_model_cached_nodes(p_num, t_num)
    return _pnjl_model_log_sum(masses, p_mesh, cosθ_mesh, coefficients, Φ, Φbar, mu_vec, T_fm, xi)
end

"""返回 (quark, antiquark) 数密度（单位 fm⁻³）。

注意：此处不再委托 legacy `Thermodynamics.calculate_number_densities`，
而是在 models 侧直接按分布函数进行数值积分（便于后续彻底解耦 legacy core）。
"""
function number_densities(
    model::PNJLModel,
    x_state,
    T_fm,
    mu_vec;
    thermal_nodes=nothing,
    p_num::Int=PNJLIntegrals.DEFAULT_MOMENTUM_COUNT,
    t_num::Int=PNJLIntegrals.DEFAULT_THETA_COUNT,
    xi=0.0,
    kwargs...
)
    st = x_state isa MeanFieldState ? x_state : MeanFieldState(x_state)
    Φ = Float64(st.Phi)
    Φbar = Float64(st.PhiBar)
    T_val = Float64(T_fm)
    xi_val = Float64(xi)

    φ = SVector{3, Float64}(Float64(st.phi[1]), Float64(st.phi[2]), Float64(st.phi[3]))
    masses = calculate_mass_vec(model, φ)

    nodes = isnothing(thermal_nodes) ? _pnjl_model_cached_nodes(p_num, t_num) : thermal_nodes
    thermal_p_mesh, cosθ_mesh, thermal_coefficients = nodes
    pref = 2 * Int(model.consts.N_color)

    acc_q = MVector{3, Float64}(0.0, 0.0, 0.0)
    acc_aq = MVector{3, Float64}(0.0, 0.0, 0.0)

    mu_vec = normalize_mu_vec(mu_vec)

    @inbounds for i in 1:3
        mass_i = Float64(masses[i])
        mu_i = Float64(mu_vec[i])
        total_q = 0.0
        total_aq = 0.0
        for idx in eachindex(thermal_p_mesh)
            p = Float64(thermal_p_mesh[idx])
            cosθ = Float64(cosθ_mesh[idx])
            w = Float64(thermal_coefficients[idx])
            total_q += w * pref * pnjl_quark_distribution_aniso(p, mass_i, mu_i, T_val, Φ, Φbar, xi_val, cosθ)
            total_aq += w * pref * pnjl_antiquark_distribution_aniso(p, mass_i, mu_i, T_val, Φ, Φbar, xi_val, cosθ)
        end
        acc_q[i] = total_q
        acc_aq[i] = total_aq
    end

    return (quark=SVector{3}(acc_q), antiquark=SVector{3}(acc_aq))
end

"""solve_gap(::PNJLModel, T, mu_vec; kwargs...) -> MeanFieldState

阶段 0：先复用 legacy 求解器（src/pnjl/solver/Solver.jl），将其结果适配到 models 的 `MeanFieldState`。

当前限制：仅支持对称化学势 `mu_u == mu_d == mu_s`（因为 legacy FixedMu 模式以标量 μ 为入口）。
"""
function solve_gap(
    model::PNJLModel,
    T_fm,
    mu_vec;
    solver_backend::Symbol=:legacy,
    xi::Real=0.0,
    p_num::Int=PNJLIntegrals.DEFAULT_MOMENTUM_COUNT,
    t_num::Int=PNJLIntegrals.DEFAULT_THETA_COUNT,
    thermo_backend::Symbol=:models,
    kwargs...
)
    if solver_backend === :legacy
        _ensure_pnjl_solver_loaded()

        μ = normalize_mu_vec(mu_vec)
        if !(μ[1] == μ[2] == μ[3])
            throw(ArgumentError("PNJLModel.solve_gap(::legacy) requires mu_u == mu_d == mu_s (FixedMu mode). Got mu_vec=$μ"))
        end

        # NOTE: We may have just `include`d the legacy solver module; use `invokelatest`
        # to avoid world-age issues when running inside tests or AD contexts.
        mode = Base.invokelatest(Main.Solver.FixedMu)
        res = Base.invokelatest(Main.Solver.solve, mode, T_fm, μ[1];
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            thermo_backend=thermo_backend,
            kwargs...)

        res.converged || error("PNJLModel.solve_gap did not converge (residual_norm=$(res.residual_norm))")
        return MeanFieldState(res.x_state)
    elseif solver_backend === :models
        # Use the generic 5D solver (gap_residual = ∇Ω). Do NOT pass legacy-only
        # kwargs (e.g. thermo_backend) down to gap_residual/omega.
        return invoke(
            solve_gap,
            Tuple{AbstractPNJLModel, Any, Any},
            model,
            T_fm,
            mu_vec;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            kwargs...,
        )
    else
        throw(ArgumentError("unknown solver_backend=$solver_backend (expected :legacy or :models)"))
    end
end
