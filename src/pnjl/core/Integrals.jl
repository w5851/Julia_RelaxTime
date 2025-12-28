"""
    Integrals

PNJL 模型积分计算模块，包含真空项和热项的积分实现。

## 公式来源
- 各向同性公式：docs/reference/formula/pnjl/Omega_各向同性.md
- 各向异性公式：docs/reference/formula/pnjl/Omega_RS各向异性.md

## 主要功能
- 积分节点缓存管理
- 真空项解析积分（式 2.26 真空部分）
- 热项数值积分（支持 RS 各向异性变形）

## 关键公式
真空项解析积分：
```math
I(Λ, M) = \\frac{1}{16π²}\\left[Λ\\sqrt{Λ² + M²}(2Λ² + M²) - M⁴\\ln\\frac{Λ + \\sqrt{Λ² + M²}}{M}\\right]
"""
module Integrals

using Base.MathConstants: π
using StaticArrays

# 加载 GaussLegendre 模块
const _GAUSSLEGENDRE_PATH = normpath(joinpath(@__DIR__, "..", "..", "integration", "GaussLegendre.jl"))
if !isdefined(Main, :GaussLegendre)
    Base.include(Main, _GAUSSLEGENDRE_PATH)
end

using Main.GaussLegendre:
    gauleg,
    DEFAULT_COSΘ_HALF_NODES,
    DEFAULT_COSΘ_HALF_WEIGHTS,
    DEFAULT_MOMENTUM_NODES,
    DEFAULT_MOMENTUM_WEIGHTS

# 常量导入 - 使用绝对路径加载
const _CONSTANTS_PATH = normpath(joinpath(@__DIR__, "..", "..", "Constants_PNJL.jl"))
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PATH)
end
using Main.Constants_PNJL: Λ_inv_fm, N_color

export cached_nodes, vacuum_integral, calculate_energy_sum, calculate_log_sum
export DEFAULT_THETA_COUNT, DEFAULT_MOMENTUM_COUNT

# ============================================================================
# 积分节点常量与缓存
# ============================================================================

const DEFAULT_THETA_COUNT = length(DEFAULT_COSΘ_HALF_NODES)
const DEFAULT_MOMENTUM_COUNT = length(DEFAULT_MOMENTUM_NODES)
const THETA_DEFAULT_NODES = DEFAULT_COSΘ_HALF_NODES
const THETA_DEFAULT_WEIGHTS = DEFAULT_COSΘ_HALF_WEIGHTS .* 2.0
const THERMAL_DEFAULT_NODES = DEFAULT_MOMENTUM_NODES
const THERMAL_DEFAULT_WEIGHTS = DEFAULT_MOMENTUM_WEIGHTS

"""积分节点缓存：(p_num, t_num) -> (p_mesh, cosθ_mesh, coefficients)"""
const NODE_CACHE = Dict{Tuple{Int, Int}, NTuple{3, Matrix{Float64}}}()

# ============================================================================
# 节点生成函数
# ============================================================================

@inline function theta_nodes(t_num::Int)
    if t_num == DEFAULT_THETA_COUNT
        return THETA_DEFAULT_NODES, THETA_DEFAULT_WEIGHTS
    end
    nodes, weights = gauleg(0.0, 1.0, t_num)
    return nodes, weights .* 2.0
end

@inline function thermal_nodes(p_num::Int)
    if p_num == DEFAULT_MOMENTUM_COUNT
        return THERMAL_DEFAULT_NODES, THERMAL_DEFAULT_WEIGHTS
    end
    return gauleg(0.0, 10.0, p_num)
end

"""
    build_nodes(p_num, t_num) -> (p_mesh, cosθ_mesh, coefficients)

构建二维积分网格节点和权重系数。
"""
function build_nodes(p_num::Int, t_num::Int)
    momentum_nodes, momentum_weights = thermal_nodes(p_num)
    cosθ_nodes, cosθ_weights = theta_nodes(t_num)

    thermal_p_mesh = repeat(momentum_nodes, 1, t_num)
    cosθ_mesh = repeat(cosθ_nodes', p_num, 1)
    weight_mesh = momentum_weights * cosθ_weights'
    thermal_coefficients = weight_mesh .* thermal_p_mesh.^2 ./ (2 * π)^2

    return (thermal_p_mesh, cosθ_mesh, thermal_coefficients)
end

"""
    cached_nodes(p_num, t_num) -> (p_mesh, cosθ_mesh, coefficients)

获取缓存的积分节点，若不存在则创建并缓存。
"""
function cached_nodes(p_num::Int, t_num::Int)
    key = (p_num, t_num)
    return get!(NODE_CACHE, key) do
        build_nodes(p_num, t_num)
    end
end

# ============================================================================
# 真空项积分
# ============================================================================

"""
    vacuum_integral(mass) -> Float64

计算单个味道的真空项解析积分：
```math
I(Λ, M) = \\frac{1}{16π²}\\left[Λ\\sqrt{Λ² + M²}(2Λ² + M²) - M⁴\\ln\\frac{Λ + \\sqrt{Λ² + M²}}{M}\\right]
```
"""
@inline function vacuum_integral(mass::TF) where {TF}
    Λ = convert(TF, Λ_inv_fm)
    mass_abs = abs(mass)
    epsilon = one(mass_abs) * 1e-12
    mass_safe = mass_abs + epsilon
    sqrt_term = sqrt(Λ^2 + mass_safe^2)
    poly_part = Λ * sqrt_term * (2 * Λ^2 + mass_safe^2)
    log_term = mass_safe^4 * log((Λ + sqrt_term) / mass_safe)
    return (poly_part - log_term) / (16 * π^2)
end

"""
    calculate_energy_sum(masses) -> Float64

计算真空能量项：-2Nc ∑_i I(Λ, M_i)
"""
function calculate_energy_sum(masses::SVector{3, TF}) where {TF}
    total = zero(TF)
    @inbounds for i in 1:3
        total += vacuum_integral(masses[i])
    end
    return -2 * N_color * total
end

# ============================================================================
# 热项积分
# ============================================================================

"""
    safe_log(x; min_val=1e-16) -> Float64

安全的对数函数，避免 log(0) 或 log(负数)。
"""
@inline function safe_log(x; min_val=1e-16)
    x <= 0 && return log(min_val)
    return x < min_val ? log(min_val) : log(x)
end

"""
    calculate_energy_isotropic(mass, p) -> Float64

各向同性色散关系：E = √(p² + m²)
"""
@inline function calculate_energy_isotropic(mass_i, p)
    return sqrt(p^2 + mass_i^2)
end

"""
    calculate_energy_anisotropic(mass, p, xi, cosθ) -> Float64

各向异性色散关系（RS 变形）：E = √(p² + m² + ξ(p·cosθ)²)
"""
@inline function calculate_energy_anisotropic(mass_i, p, xi, t)
    return sqrt(p^2 + mass_i^2 + xi * (p * t)^2)
end

"""
    calculate_log_term(E, mu, T, Φ, Φ̄) -> Float64

计算 Polyakov loop 修正的对数项：ln(f₊) + ln(f₋)
"""
@inline function calculate_log_term(E_i, mu_i, T_fm, Φ, Φ̄)
    invT = 1.0 / T_fm
    x_i = (E_i - mu_i) * invT
    x_i_anti = (E_i + mu_i) * invT

    exp1 = exp(-x_i)
    exp2 = exp1 * exp1
    exp3 = exp1 * exp2
    exp1_anti = exp(-x_i_anti)
    exp2_anti = exp1_anti * exp1_anti
    exp3_anti = exp1_anti * exp2_anti

    f1_val = 1.0 + 3.0 * Φ * exp1 + 3.0 * Φ̄ * exp2 + exp3
    f2_val = 1.0 + 3.0 * Φ̄ * exp1_anti + 3.0 * Φ * exp2_anti + exp3_anti

    return safe_log(f1_val) + safe_log(f2_val)
end

"""
    calculate_log_sum(masses, p_nodes, cosθ_nodes, coefficients, Φ, Φ̄, mu_vec, T, xi) -> Float64

计算热项对数和（各向异性）：
```math
-2T \\sum_i \\int \\frac{d³p}{(2π)³} [\\mathcal{Q}_1^{aniso} + \\mathcal{Q}_2^{aniso}]
```
"""
function calculate_log_sum(masses::SVector{3, TF}, p_nodes, cosθ_nodes, coefficients, Φ, Φ̄, mu_vec, T_fm, xi) where {TF}
    total = zero(TF)
    @inbounds for i in 1:3
        mass_i = masses[i]
        mu_i = mu_vec[i]
        for idx in eachindex(p_nodes)
            E_i = calculate_energy_anisotropic(mass_i, p_nodes[idx], xi, cosθ_nodes[idx])
            total += calculate_log_term(E_i, mu_i, T_fm, Φ, Φ̄) * coefficients[idx]
        end
    end
    return -2 * T_fm * total
end

end # module Integrals
