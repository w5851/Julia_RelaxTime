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

# 导入 ConstraintModes
include("ConstraintModes.jl")
using .ConstraintModes: ConstraintMode, FixedMu, FixedRho, FixedEntropy, FixedSigma, state_dim

# 导入 core 模块
include(joinpath(@__DIR__, "..", "core", "Thermodynamics.jl"))
using .Thermodynamics: calculate_pressure, calculate_rho, calculate_thermo, ρ0
using .Thermodynamics.Integrals: cached_nodes

export gap_conditions, build_conditions, build_residual!
export GapParams

# ============================================================================
# 参数结构
# ============================================================================

"""
    GapParams

能隙方程求解所需的参数集合。

# 字段
- `T_fm::Float64`: 温度 (fm⁻¹)
- `thermal_nodes`: 积分节点
- `xi::Float64`: 各向异性参数
"""
struct GapParams{TN}
    T_fm::Float64
    thermal_nodes::TN
    xi::Float64
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
        calculate_pressure(y_s, mu_vec, params.T_fm, params.thermal_nodes, params.xi)
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
        local_params = GapParams(T_fm, params.thermal_nodes, params.xi)
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
        local_params = GapParams(T_fm, params.thermal_nodes, params.xi)
        
        # 5 个能隙方程
        gap = gap_conditions(x_state, mu_vec, local_params)
        
        # 化学势相等约束
        mu_eq1 = x[6] - x[7]  # μ_u = μ_d
        mu_eq2 = x[7] - x[8]  # μ_d = μ_s
        
        # 密度约束
        rho_vec = calculate_rho(x_state, mu_vec, T_fm, params.thermal_nodes, params.xi)
        rho_constraint = sum(rho_vec) / (3.0 * ρ0) - mode.rho_target
        
        return [gap..., mu_eq1, mu_eq2, rho_constraint]
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
        local_params = GapParams(T_fm, params.thermal_nodes, params.xi)
        
        # 5 个能隙方程
        gap = gap_conditions(x_state, mu_vec, local_params)
        
        # 化学势相等约束
        mu_eq1 = x[6] - x[7]
        mu_eq2 = x[7] - x[8]
        
        # 熵密度约束
        _, _, entropy, _ = calculate_thermo(x_state, mu_vec, T_fm, params.thermal_nodes, params.xi)
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
        local_params = GapParams(T_fm, params.thermal_nodes, params.xi)
        
        # 5 个能隙方程
        gap = gap_conditions(x_state, mu_vec, local_params)
        
        # 化学势相等约束
        mu_eq1 = x[6] - x[7]
        mu_eq2 = x[7] - x[8]
        
        # 比熵约束 σ = s/n_B
        _, rho_norm, entropy, _ = calculate_thermo(x_state, mu_vec, T_fm, params.thermal_nodes, params.xi)
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
        rho = calculate_rho(x_state, mu_state, params.T_fm, params.thermal_nodes, params.xi)
        F[8] = sum(rho) / (3.0 * ρ0) - mode.rho_target
        
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
        _, _, entropy, _ = calculate_thermo(x_state, mu_state, params.T_fm, params.thermal_nodes, params.xi)
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
        _, rho_norm, entropy, _ = calculate_thermo(x_state, mu_state, params.T_fm, params.thermal_nodes, params.xi)
        n_B = rho_norm * ρ0
        sigma = n_B > 1e-12 ? entropy / n_B : 0.0
        F[8] = sigma - mode.sigma_target
        
        return nothing
    end
end

end # module Conditions
