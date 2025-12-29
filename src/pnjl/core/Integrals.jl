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
export calculate_log_sum_derivatives
export DEFAULT_THETA_COUNT, DEFAULT_MOMENTUM_COUNT
export calculate_energy_isotropic, calculate_energy_anisotropic

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

# ============================================================================
# 热项导数的解析计算（避免嵌套 ForwardDiff）
# ============================================================================

"""
    calculate_log_term_derivatives(E, mu, T, Φ, Φ̄) -> (log_term, d_log_dmu, d_log_dT)

计算 Polyakov loop 修正的对数项及其对 μ 和 T 的导数。

返回：
- `log_term`: ln(f₊) + ln(f₋)
- `d_log_dmu`: ∂[ln(f₊) + ln(f₋)]/∂μ
- `d_log_dT`: ∂[ln(f₊) + ln(f₋)]/∂T
"""
@inline function calculate_log_term_derivatives(E_i, mu_i, T_fm, Φ, Φ̄)
    invT = 1.0 / T_fm
    x = (E_i - mu_i) * invT  # 夸克
    y = (E_i + mu_i) * invT  # 反夸克

    # 指数项
    exp1_x = exp(-x)
    exp2_x = exp1_x * exp1_x
    exp3_x = exp1_x * exp2_x
    exp1_y = exp(-y)
    exp2_y = exp1_y * exp1_y
    exp3_y = exp1_y * exp2_y

    # f₊ 和 f₋
    f_plus = 1.0 + 3.0 * Φ * exp1_x + 3.0 * Φ̄ * exp2_x + exp3_x
    f_minus = 1.0 + 3.0 * Φ̄ * exp1_y + 3.0 * Φ * exp2_y + exp3_y

    # 对数项
    log_term = safe_log(f_plus) + safe_log(f_minus)

    # ∂f₊/∂μ = (1/T) * (3Φ·e^{-x} + 6Φ̄·e^{-2x} + 3e^{-3x})
    # ∂f₋/∂μ = (-1/T) * (3Φ̄·e^{-y} + 6Φ·e^{-2y} + 3e^{-3y})
    df_plus_dmu = invT * (3.0 * Φ * exp1_x + 6.0 * Φ̄ * exp2_x + 3.0 * exp3_x)
    df_minus_dmu = -invT * (3.0 * Φ̄ * exp1_y + 6.0 * Φ * exp2_y + 3.0 * exp3_y)

    # ∂ln(f)/∂μ = (1/f) * ∂f/∂μ
    d_log_dmu = df_plus_dmu / f_plus + df_minus_dmu / f_minus

    # ∂f₊/∂T = (x/T) * (3Φ·e^{-x} + 6Φ̄·e^{-2x} + 3e^{-3x})
    # ∂f₋/∂T = (y/T) * (3Φ̄·e^{-y} + 6Φ·e^{-2y} + 3e^{-3y})
    df_plus_dT = (x * invT) * (3.0 * Φ * exp1_x + 6.0 * Φ̄ * exp2_x + 3.0 * exp3_x)
    df_minus_dT = (y * invT) * (3.0 * Φ̄ * exp1_y + 6.0 * Φ * exp2_y + 3.0 * exp3_y)

    # ∂ln(f)/∂T = (1/f) * ∂f/∂T
    d_log_dT = df_plus_dT / f_plus + df_minus_dT / f_minus

    return (log_term, d_log_dmu, d_log_dT)
end

"""
    calculate_log_sum_derivatives(masses, p_nodes, cosθ_nodes, coefficients, Φ, Φ̄, mu_vec, T, xi)
        -> (log_sum, d_log_sum_dmu::SVector{3}, d_log_sum_dT)

计算热项对数和及其对 μ_i 和 T 的导数（不使用 ForwardDiff）。

log_sum = -2T ∑_i ∫ [ln(f₊) + ln(f₋)] d³p/(2π)³

返回：
- `log_sum`: 热项对数和
- `d_log_sum_dmu`: ∂log_sum/∂μ_i (i=1,2,3)
- `d_log_sum_dT`: ∂log_sum/∂T
"""
function calculate_log_sum_derivatives(masses::SVector{3, TF}, p_nodes, cosθ_nodes, coefficients, 
                                       Φ, Φ̄, mu_vec, T_fm, xi) where {TF}
    # 确定输出类型（支持 Dual 类型）
    RT = promote_type(TF, typeof(T_fm), eltype(mu_vec))
    
    # 累加器
    total_log = zero(RT)
    total_dmu_1 = zero(RT)
    total_dmu_2 = zero(RT)
    total_dmu_3 = zero(RT)
    total_dT = zero(RT)

    @inbounds for i in 1:3
        mass_i = masses[i]
        mu_i = mu_vec[i]
        flavor_log = zero(RT)
        flavor_dmu = zero(RT)
        flavor_dT = zero(RT)

        for idx in eachindex(p_nodes)
            p = p_nodes[idx]
            cosθ = cosθ_nodes[idx]
            w = coefficients[idx]
            E_i = calculate_energy_anisotropic(mass_i, p, xi, cosθ)
            
            log_term, d_log_dmu, d_log_dT = calculate_log_term_derivatives(E_i, mu_i, T_fm, Φ, Φ̄)
            
            flavor_log += w * log_term
            flavor_dmu += w * d_log_dmu
            flavor_dT += w * d_log_dT
        end

        total_log += flavor_log
        if i == 1
            total_dmu_1 = flavor_dmu
        elseif i == 2
            total_dmu_2 = flavor_dmu
        else
            total_dmu_3 = flavor_dmu
        end
        total_dT += flavor_dT
    end

    # log_sum = -2T * total_log
    # ∂log_sum/∂μ_i = -2T * ∂total_log/∂μ_i
    # ∂log_sum/∂T = -2 * total_log - 2T * ∂total_log/∂T
    log_sum = -2 * T_fm * total_log
    d_log_sum_dmu = SVector{3, RT}(-2 * T_fm * total_dmu_1, 
                                   -2 * T_fm * total_dmu_2, 
                                   -2 * T_fm * total_dmu_3)
    d_log_sum_dT = -2 * total_log - 2 * T_fm * total_dT

    return (log_sum, d_log_sum_dmu, d_log_sum_dT)
end

end # module Integrals
