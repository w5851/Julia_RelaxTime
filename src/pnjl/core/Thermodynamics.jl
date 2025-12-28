"""
    Thermodynamics

PNJL 模型热力学量计算模块。

## 主要功能
- 夸克有效质量计算
- 手征凝聚项计算
- Polyakov loop 有效势计算
- 巨热力学势 Ω 计算
- 压强、熵密度、能量密度、粒子数密度计算

## 单位约定
所有量使用自然单位 (fm⁻¹)，温度和化学势需要从 MeV 转换。
"""
module Thermodynamics

using Base.MathConstants: π
using StaticArrays
using ForwardDiff

# 常量导入 - 使用绝对路径加载
const _CONSTANTS_PATH = normpath(joinpath(@__DIR__, "..", "..", "Constants_PNJL.jl"))
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PATH)
end
using Main.Constants_PNJL:
    ħc_MeV_fm,
    N_color,
    ρ0_inv_fm3,
    m_ud0_inv_fm,
    m_s0_inv_fm,
    G_fm2,
    K_fm5,
    T0_inv_fm,
    a0,
    a1,
    a2,
    b3

# 导入 Integrals 模块
include("Integrals.jl")
using .Integrals:
    safe_log,
    calculate_energy_sum,
    calculate_log_sum,
    calculate_energy_isotropic,
    calculate_energy_anisotropic

export calculate_mass_vec, calculate_chiral, calculate_U
export calculate_pressure, calculate_omega
export calculate_rho, calculate_thermo, calculate_number_densities
export ρ0

# 导出常用常量
const ρ0 = ρ0_inv_fm3

# ============================================================================
# 质量计算
# ============================================================================

"""
    calculate_mass_vec(φ) -> SVector{3}

计算三味夸克的有效质量：
```math
M_i = m_{i0} - 4Gφ_i + 2Kφ_jφ_k
```

# 参数
- `φ`: 手征凝聚 [φ_u, φ_d, φ_s]

# 返回
- SVector{3}: 有效质量 [M_u, M_d, M_s] (fm⁻¹)
"""
@inline function calculate_mass_vec(φ::SVector{3, TF}) where {TF}
    φ_u, φ_d, φ_s = φ
    return SVector{3, TF}(
        m_ud0_inv_fm - 4 * G_fm2 * φ_u + 2 * K_fm5 * φ_d * φ_s,
        m_ud0_inv_fm - 4 * G_fm2 * φ_d + 2 * K_fm5 * φ_u * φ_s,
        m_s0_inv_fm - 4 * G_fm2 * φ_s + 2 * K_fm5 * φ_u * φ_d,
    )
end

"""
    calculate_mass_vec(x_state) -> SVector{3}

从完整状态向量提取 φ 并计算质量。
"""
@inline function calculate_mass_vec(x_state::SVector{5, TF}) where {TF}
    φ = SVector{3, TF}(x_state[1], x_state[2], x_state[3])
    return calculate_mass_vec(φ)
end

# ============================================================================
# 势能项计算
# ============================================================================

"""
    calculate_chiral(φ) -> Float64

计算手征凝聚项：
```math
χ = 2G \\sum_i φ_i² - 4K φ_u φ_d φ_s
```
"""
@inline function calculate_chiral(φ::SVector{3, TF}) where {TF}
    2 * G_fm2 * sum(φ .^ 2) - 4 * K_fm5 * prod(φ)
end

"""
    calculate_U(T, Φ, Φ̄) -> Float64

计算 Polyakov loop 有效势（对数形式）：
```math
\\frac{U}{T⁴} = -\\frac{a(T)}{2}Φ\\bar{Φ} + b(T)\\ln[1 - 6Φ\\bar{Φ} + 4(Φ³+\\bar{Φ}³) - 3(Φ\\bar{Φ})²]
```
"""
@inline function calculate_U(T_fm::Real, Φ::Real, Φ̄::Real)
    T_ratio = T0_inv_fm / T_fm
    Ta = a0 + a1 * T_ratio + a2 * T_ratio^2
    Tb = b3 * T_ratio^3
    value = 1 - 6 * Φ̄ * Φ + 4 * (Φ̄^3 + Φ^3) - 3 * (Φ̄ * Φ)^2
    return T_fm^4 * (-0.5 * Ta * Φ̄ * Φ + Tb * safe_log(value))
end

# ============================================================================
# 巨热力学势与压强
# ============================================================================

"""
    calculate_omega(x_state, mu_vec, T, thermal_nodes, xi) -> Float64

计算巨热力学势 Ω：
```math
Ω = χ + U - 2N_c \\sum_i I(Λ, M_i) - 2T \\sum_i \\int [\\mathcal{Q}_1 + \\mathcal{Q}_2]
```

# 参数
- `x_state`: 状态向量 [φ_u, φ_d, φ_s, Φ, Φ̄]
- `mu_vec`: 化学势向量 [μ_u, μ_d, μ_s] (fm⁻¹)
- `T`: 温度 (fm⁻¹)
- `thermal_nodes`: 积分节点 (p_mesh, cosθ_mesh, coefficients)
- `xi`: 各向异性参数
"""
function calculate_omega(x_state::SVector{5, TF}, mu_vec::AbstractVector{TM}, T_fm::TR, thermal_nodes, xi) where {TF, TM, TR}
    φ = SVector{3, TF}(x_state[1], x_state[2], x_state[3])
    Φ, Φ̄ = x_state[4], x_state[5]

    chi = calculate_chiral(φ)
    U = calculate_U(T_fm, Φ, Φ̄)
    masses = calculate_mass_vec(φ)

    thermal_p_mesh, cosθ_mesh, thermal_coefficients = thermal_nodes

    energy_sum = calculate_energy_sum(masses)
    log_sum = calculate_log_sum(masses, thermal_p_mesh, cosθ_mesh, thermal_coefficients, Φ, Φ̄, mu_vec, T_fm, xi)
    
    return chi + U + energy_sum + log_sum
end

"""
    calculate_pressure(x_state, mu_vec, T, thermal_nodes, xi) -> Float64

计算压强：P = -Ω
"""
function calculate_pressure(x_state::SVector{5, TF}, mu_vec::AbstractVector{TM}, T_fm::TR, thermal_nodes, xi) where {TF, TM, TR}
    return -calculate_omega(x_state, mu_vec, T_fm, thermal_nodes, xi)
end

# ============================================================================
# 热力学量计算
# ============================================================================

"""
    calculate_rho(x_state, mu_vec, T, thermal_nodes, xi) -> SVector{3}

计算粒子数密度：ρ_i = -∂Ω/∂μ_i = ∂P/∂μ_i
"""
function calculate_rho(x_state::SVector{5, TF}, mu_vec::AbstractVector{TM}, T_fm::TR, thermal_nodes, xi) where {TF, TM, TR}
    pressure_mu = μ -> calculate_pressure(x_state, μ, T_fm, thermal_nodes, xi)
    grad = ForwardDiff.gradient(pressure_mu, mu_vec)
    grad_type = typeof(grad[1])
    return SVector{3, grad_type}(Tuple(grad))
end

"""
    calculate_thermo(x_state, mu_vec, T, thermal_nodes, xi) -> (P, ρ_norm, s, ε)

计算所有热力学量。

# 返回
- `P`: 压强 (fm⁻⁴)
- `ρ_norm`: 归一化重子数密度 ρ_B/ρ₀
- `s`: 熵密度 (fm⁻³)
- `ε`: 能量密度 (fm⁻⁴)
"""
function calculate_thermo(x_state::SVector{5, TF}, mu_vec::AbstractVector{TM}, T_fm::TR, thermal_nodes, xi) where {TF, TM, TR}
    rho_vec = calculate_rho(x_state, mu_vec, T_fm, thermal_nodes, xi)
    rho_norm = sum(rho_vec) / (3.0 * ρ0)
    
    pressure_T = τ -> calculate_pressure(x_state, mu_vec, τ, thermal_nodes, xi)
    entropy = ForwardDiff.derivative(pressure_T, T_fm)
    
    pressure = calculate_pressure(x_state, mu_vec, T_fm, thermal_nodes, xi)
    energy = -pressure + sum(mu_vec .* rho_vec) + T_fm * entropy
    
    return pressure, rho_norm, entropy, energy
end

"""
    calculate_number_densities(x_state, mu_vec, T, thermal_nodes, xi) -> (quark, antiquark)

分别计算夸克与反夸克的数密度（各向异性分布）。

# 返回
- `quark`: 夸克数密度 [n_u, n_d, n_s]
- `antiquark`: 反夸克数密度 [n_ū, n_d̄, n_s̄]
"""
function calculate_number_densities(x_state::SVector{5, TF}, mu_vec::AbstractVector{TM}, T_fm::TR, thermal_nodes, xi) where {TF, TM, TR}
    # 需要导入分布函数
    include_path = joinpath(@__DIR__, "..", "..", "QuarkDistribution_Aniso.jl")
    if !isdefined(@__MODULE__, :PNJLQuarkDistributions_Aniso)
        include(include_path)
    end
    
    φ = SVector{3, TF}(x_state[1], x_state[2], x_state[3])
    Φ, Φ̄ = x_state[4], x_state[5]
    masses = calculate_mass_vec(φ)

    thermal_p_mesh, cosθ_mesh, thermal_coefficients = thermal_nodes
    pref = 2 * N_color
    acc_q = MVector{3, promote_type(TF, TM, eltype(thermal_p_mesh))}(0, 0, 0)
    acc_aq = similar(acc_q)

    @inbounds for i in 1:3
        mass_i = masses[i]
        mu_i = mu_vec[i]
        total_q = zero(acc_q[i])
        total_aq = zero(acc_aq[i])
        for idx in eachindex(thermal_p_mesh)
            p = thermal_p_mesh[idx]
            cosθ = cosθ_mesh[idx]
            w = thermal_coefficients[idx]
            total_q += w * pref * Main.PNJLQuarkDistributions_Aniso.quark_distribution_aniso(p, mass_i, mu_i, T_fm, Φ, Φ̄, xi, cosθ)
            total_aq += w * pref * Main.PNJLQuarkDistributions_Aniso.antiquark_distribution_aniso(p, mass_i, mu_i, T_fm, Φ, Φ̄, xi, cosθ)
        end
        acc_q[i] = total_q
        acc_aq[i] = total_aq
    end
    return (quark = SVector{3}(acc_q), antiquark = SVector{3}(acc_aq))
end

end # module Thermodynamics
