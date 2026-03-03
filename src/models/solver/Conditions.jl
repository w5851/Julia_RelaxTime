"""
    Conditions

PNJL 求解器条件函数定义。

## 核心思想
将能隙方程的条件定义与求解器分离，支持多种约束模式。

## 主要函数
- `gap_conditions`: 核心能隙方程（5维，所有模式共用）
- `build_conditions`: 根据模式构建完整条件函数
- `build_residual!`: 构建 NLsolve 兼容的残差函数
"""
module Conditions

using StaticArrays
using ForwardDiff

# 从 Models 域导入
import Main.Models: ConstraintMode, FixedMu, FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma, state_dim

const _PNJL_CORE_PATH = normpath(joinpath(@__DIR__, "..", "pnjl_physics", "PNJLCore.jl"))
if !isdefined(@__MODULE__, :PNJLCore)
    include(_PNJL_CORE_PATH)
end
using .PNJLCore: cached_nodes

# Unified thermo facade (legacy vs models)
const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

const _CONSTANTS_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "Constants_PNJL.jl"))
IncludeOnce.include_once!(Main, :Constants_PNJL, _CONSTANTS_PATH)
using Main.Constants_PNJL: ρ0_inv_fm3
const ρ0 = ρ0_inv_fm3

const _THERMO_FACADE_PATH = normpath(joinpath(@__DIR__, "..", "pnjl_physics", "core", "ThermoFacade.jl"))
const ThermoFacade = IncludeOnce.include_once!(Main, :ThermoFacade, _THERMO_FACADE_PATH)

import Main.Models: AbstractQCDModel, AbstractPNJLModel, PNJLModel, PNJLMagneticModel, RPNJLModel

export gap_conditions, build_conditions, build_residual!
export GapParams

@inline _model_kind_symbol(::AbstractPNJLModel) = :PNJL
@inline _model_kind_symbol(::PNJLModel) = :PNJL
@inline _model_kind_symbol(::PNJLMagneticModel) = :PNJL
@inline _model_kind_symbol(::RPNJLModel) = :RPNJL
@inline _model_kind_symbol(::AbstractQCDModel) = :PNJL

# ============================================================================
# 参数结构
# ============================================================================

"""
    GapParams

能隙方程求解所需的参数集合。

# 字段
- `T_fm`: 温度 (fm⁻¹)
- `thermal_nodes`: 积分节点
- `xi`: 各向异性参数
"""
struct GapParams{TT, TN, TX}
    T_fm::TT
    thermal_nodes::TN
    xi::TX

    p_num::Int
    t_num::Int
    model_kind::Symbol
end

@inline function GapParams(
    T_fm::TT,
    thermal_nodes::TN,
    xi::TX;
    p_num::Int=size(thermal_nodes[1], 1),
    t_num::Int=size(thermal_nodes[1], 2),
    model_kind::Symbol=:PNJL,
) where {TT, TN, TX}
    return GapParams(T_fm, thermal_nodes, xi, p_num, t_num, model_kind)
end

@inline function _rho_vec(x_state, mu_vec, T_fm, params::GapParams)
    return ThermoFacade.calculate_rho_backend(
        x_state,
        mu_vec,
        T_fm;
        model_kind=params.model_kind,
        p_num=params.p_num,
        t_num=params.t_num,
        thermal_nodes=params.thermal_nodes,
        xi=params.xi,
    )
end

@inline function _thermo_tuple(x_state, mu_vec, T_fm, params::GapParams)
    return ThermoFacade.calculate_thermo_backend(
        x_state,
        mu_vec,
        T_fm;
        model_kind=params.model_kind,
        p_num=params.p_num,
        t_num=params.t_num,
        thermal_nodes=params.thermal_nodes,
        xi=params.xi,
    )
end

@inline function _safe_density_ratio(num, den; eps::Float64=1e-12)
    abs(den) >= eps && return num / den
    den_safe = den >= 0 ? eps : -eps
    return num / den_safe
end

# ============================================================================
# 核心能隙条件（所有模式共用）
# ============================================================================

"""
    gap_conditions(x_state, mu_vec, params) -> SVector{5}

计算 5 个能隙方程的残差：∂Ω/∂φ_u, ∂Ω/∂φ_d, ∂Ω/∂φ_s, ∂Ω/∂Φ, ∂Ω/∂Φ̄

这是所有求解模式的核心，通过对压强求梯度实现。

# 参数
- `x_state`: 状态向量 [φ_u, φ_d, φ_s, Φ, Φ̄]
- `mu_vec`: 化学势向量 [μ_u, μ_d, μ_s]
- `params`: GapParams 参数结构

# 返回
- SVector{5}: 能隙方程残差（平衡态时为零）
"""
function gap_conditions(x_state::SVector{5, TF}, mu_vec::AbstractVector{TM}, params::GapParams) where {TF, TM}
    pressure_fn = y -> begin
        eltp = typeof(y[1])
        y_s = SVector{5, eltp}(Tuple(y))
        ThermoFacade.calculate_pressure_backend(
            y_s,
            mu_vec,
            params.T_fm;
            model_kind=params.model_kind,
            p_num=params.p_num,
            t_num=params.t_num,
            thermal_nodes=params.thermal_nodes,
            xi=params.xi,
        )
    end
    grad = ForwardDiff.gradient(pressure_fn, x_state)
    grad_type = typeof(grad[1])
    return SVector{5, grad_type}(Tuple(grad))
end

# ============================================================================
# 模式特定条件构建
# ============================================================================

"""
    build_conditions(mode::FixedMu, params::GapParams) -> Function

构建固定化学势模式的条件函数。

返回函数签名：(θ, x) -> residual
- θ = [T, μ]（参数）
- x = [φ_u, φ_d, φ_s, Φ, Φ̄]（状态变量）
"""
function build_conditions(::FixedMu, params::GapParams)
    return (θ, x) -> begin
        T_fm = θ[1]
        μ_fm = θ[2]
        mu_vec = SVector{3}(μ_fm, μ_fm, μ_fm)
        x_state = SVector{5}(Tuple(x))
        local_params = GapParams(T_fm, params.thermal_nodes, params.xi,
            p_num=params.p_num, t_num=params.t_num, model_kind=params.model_kind)
        return Vector(gap_conditions(x_state, mu_vec, local_params))
    end
end

"""
    build_conditions(mode::FixedRho, params::GapParams) -> Function

构建固定密度模式的条件函数。

返回函数签名：(θ, x) -> residual
- θ = [T]（参数）
- x = [φ_u, φ_d, φ_s, Φ, Φ̄, μ_u, μ_d, μ_s]（状态变量）
"""
function build_conditions(mode::FixedRho, params::GapParams)
    return (θ, x) -> begin
        T_fm = θ[1]
        x_state = SVector{5}(x[1], x[2], x[3], x[4], x[5])
        mu_vec = SVector{3}(x[6], x[7], x[8])
        local_params = GapParams(T_fm, params.thermal_nodes, params.xi,
            p_num=params.p_num, t_num=params.t_num, model_kind=params.model_kind)
        
        # 5 个能隙方程
        gap = gap_conditions(x_state, mu_vec, local_params)
        
        # 化学势相等约束
        mu_eq1 = x[6] - x[7]  # μ_u = μ_d
        mu_eq2 = x[7] - x[8]  # μ_d = μ_s
        
        # 密度约束
        rho_vec = _rho_vec(x_state, mu_vec, T_fm, local_params)
        rho_constraint = sum(rho_vec) / (3.0 * ρ0) - mode.rho_target
        
        return [gap..., mu_eq1, mu_eq2, rho_constraint]
    end
end

"""
    build_conditions(mode::FixedAsymmetricRho, params::GapParams) -> Function

构建固定非对称约束模式的条件函数。

返回函数签名：(θ, x) -> residual
- θ = [T]（参数）
- x = [φ_u, φ_d, φ_s, Φ, Φ̄, μ_u, μ_d, μ_s]（状态变量）
"""
function build_conditions(mode::FixedAsymmetricRho, params::GapParams)
    return (θ, x) -> begin
        T_fm = θ[1]
        x_state = SVector{5}(x[1], x[2], x[3], x[4], x[5])
        mu_vec = SVector{3}(x[6], x[7], x[8])
        local_params = GapParams(T_fm, params.thermal_nodes, params.xi,
            p_num=params.p_num, t_num=params.t_num, model_kind=params.model_kind)

        gap = gap_conditions(x_state, mu_vec, local_params)

        rho_vec = _rho_vec(x_state, mu_vec, T_fm, local_params)
        rho_u, rho_d, rho_s = rho_vec[1], rho_vec[2], rho_vec[3]

        nB_constraint = sum(rho_vec) / (3.0 * ρ0) - mode.rho_target
        ud_ratio_constraint = _safe_density_ratio(rho_u, rho_d) - mode.ud_ratio_target
        s_constraint = rho_s - mode.s_target

        return [gap..., nB_constraint, ud_ratio_constraint, s_constraint]
    end
end

"""
    build_conditions(mode::FixedEntropy, params::GapParams) -> Function

构建固定熵密度模式的条件函数。
"""
function build_conditions(mode::FixedEntropy, params::GapParams)
    return (θ, x) -> begin
        T_fm = θ[1]
        x_state = SVector{5}(x[1], x[2], x[3], x[4], x[5])
        mu_vec = SVector{3}(x[6], x[7], x[8])
        local_params = GapParams(T_fm, params.thermal_nodes, params.xi,
            p_num=params.p_num, t_num=params.t_num, model_kind=params.model_kind)
        
        # 5 个能隙方程
        gap = gap_conditions(x_state, mu_vec, local_params)
        
        # 化学势相等约束
        mu_eq1 = x[6] - x[7]
        mu_eq2 = x[7] - x[8]
        
        # 熵密度约束
        _, _, entropy, _ = _thermo_tuple(x_state, mu_vec, T_fm, local_params)
        s_constraint = entropy - mode.s_target
        
        return [gap..., mu_eq1, mu_eq2, s_constraint]
    end
end

"""
    build_conditions(mode::FixedSigma, params::GapParams) -> Function

构建固定比熵模式的条件函数。
"""
function build_conditions(mode::FixedSigma, params::GapParams)
    return (θ, x) -> begin
        T_fm = θ[1]
        x_state = SVector{5}(x[1], x[2], x[3], x[4], x[5])
        mu_vec = SVector{3}(x[6], x[7], x[8])
        local_params = GapParams(T_fm, params.thermal_nodes, params.xi,
            p_num=params.p_num, t_num=params.t_num, model_kind=params.model_kind)
        
        # 5 个能隙方程
        gap = gap_conditions(x_state, mu_vec, local_params)
        
        # 化学势相等约束
        mu_eq1 = x[6] - x[7]
        mu_eq2 = x[7] - x[8]
        
        # 比熵约束 σ = s/n_B
        _, rho_norm, entropy, _ = _thermo_tuple(x_state, mu_vec, T_fm, local_params)
        n_B = rho_norm * ρ0  # 重子数密度
        sigma = n_B > 1e-12 ? entropy / n_B : 0.0
        sigma_constraint = sigma - mode.sigma_target
        
        return [gap..., mu_eq1, mu_eq2, sigma_constraint]
    end
end

# ============================================================================
# NLsolve 兼容的残差函数构建
# ============================================================================

"""
    build_residual!(mode::FixedMu, mu_vec, params::GapParams) -> Function

构建 NLsolve 兼容的残差函数（固定化学势模式）。

返回函数签名：(F, x) -> nothing（原地修改 F）
"""
function build_residual!(::FixedMu, mu_vec::SVector{3}, params::GapParams)
    return (F, x) -> begin
        eltp = typeof(x[1])
        x_state = SVector{5, eltp}(Tuple(x))
        core_grad = gap_conditions(x_state, mu_vec, params)
        F .= core_grad
        return nothing
    end
end

"""
    build_residual!(mode::FixedRho, params::GapParams) -> Function

构建 NLsolve 兼容的残差函数（固定密度模式）。
"""
function build_residual!(mode::FixedRho, params::GapParams)
    return (F, x) -> begin
        eltp = typeof(x[1])
        x_state = SVector{5, eltp}(Tuple(x[1:5]))
        mu_state = SVector{3, eltp}(x[6], x[7], x[8])
        
        # 能隙方程
        F[1:5] = gap_conditions(x_state, mu_state, params)
        
        # 化学势相等
        F[6] = x[6] - x[7]
        F[7] = x[7] - x[8]
        
        # 密度约束
        rho = _rho_vec(x_state, mu_state, params.T_fm, params)
        F[8] = sum(rho) / (3.0 * ρ0) - mode.rho_target
        
        return nothing
    end
end

"""
    build_residual!(mode::FixedAsymmetricRho, params::GapParams) -> Function

构建 NLsolve 兼容的残差函数（固定非对称约束模式）。
"""
function build_residual!(mode::FixedAsymmetricRho, params::GapParams)
    return (F, x) -> begin
        eltp = typeof(x[1])
        x_state = SVector{5, eltp}(Tuple(x[1:5]))
        mu_state = SVector{3, eltp}(x[6], x[7], x[8])

        F[1:5] = gap_conditions(x_state, mu_state, params)

        rho = _rho_vec(x_state, mu_state, params.T_fm, params)
        rho_u, rho_d, rho_s = rho[1], rho[2], rho[3]

        F[6] = sum(rho) / (3.0 * ρ0) - mode.rho_target
        F[7] = _safe_density_ratio(rho_u, rho_d) - mode.ud_ratio_target
        F[8] = rho_s - mode.s_target

        return nothing
    end
end

"""
    build_residual!(mode::FixedEntropy, params::GapParams) -> Function

构建 NLsolve 兼容的残差函数（固定熵密度模式）。
"""
function build_residual!(mode::FixedEntropy, params::GapParams)
    return (F, x) -> begin
        eltp = typeof(x[1])
        x_state = SVector{5, eltp}(Tuple(x[1:5]))
        mu_state = SVector{3, eltp}(x[6], x[7], x[8])
        
        # 能隙方程
        F[1:5] = gap_conditions(x_state, mu_state, params)
        
        # 化学势相等
        F[6] = x[6] - x[7]
        F[7] = x[7] - x[8]
        
        # 熵密度约束
        _, _, entropy, _ = _thermo_tuple(x_state, mu_state, params.T_fm, params)
        F[8] = entropy - mode.s_target
        
        return nothing
    end
end

"""
    build_residual!(mode::FixedSigma, params::GapParams) -> Function

构建 NLsolve 兼容的残差函数（固定比熵模式）。
"""
function build_residual!(mode::FixedSigma, params::GapParams)
    return (F, x) -> begin
        eltp = typeof(x[1])
        x_state = SVector{5, eltp}(Tuple(x[1:5]))
        mu_state = SVector{3, eltp}(x[6], x[7], x[8])
        
        # 能隙方程
        F[1:5] = gap_conditions(x_state, mu_state, params)
        
        # 化学势相等
        F[6] = x[6] - x[7]
        F[7] = x[7] - x[8]
        
        # 比熵约束
        _, rho_norm, entropy, _ = _thermo_tuple(x_state, mu_state, params.T_fm, params)
        n_B = rho_norm * ρ0
        sigma = n_B > 1e-12 ? entropy / n_B : 0.0
        F[8] = sigma - mode.sigma_target
        
        return nothing
    end
end

function build_residual!(model::AbstractQCDModel, mode::FixedMu, T_fm::Real, μ_fm::Real;
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
)
    thermal_nodes = cached_nodes(p_num, t_num)
    params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi);
        p_num=p_num,
        t_num=t_num,
        model_kind=_model_kind_symbol(model),
    )
    mu_vec = SVector{3}(Float64(μ_fm), Float64(μ_fm), Float64(μ_fm))
    return build_residual!(mode, mu_vec, params)
end

function build_residual!(model::AbstractQCDModel, mode::Union{FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma}, T_fm::Real;
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
)
    thermal_nodes = cached_nodes(p_num, t_num)
    params = GapParams(Float64(T_fm), thermal_nodes, Float64(xi);
        p_num=p_num,
        t_num=t_num,
        model_kind=_model_kind_symbol(model),
    )
    return build_residual!(mode, params)
end

end # module Conditions

