module DifferentialCrossSection

"""
# DifferentialCrossSection Module

微分散射截面计算模块，连接散射矩阵元与总散射截面/弛豫时间。

## 物理公式

在相对论性 Boltzmann 动力学理论中，微分散射截面为：

dσ/dt = [1/(16π s₁₂⁺ s₁₂⁻)] · [1/(4Nc²) Σ|M|²]

其中：
- s₁₂⁺ = s - (m₁ + m₂)² （入射粒子质量和的阈值）
- s₁₂⁻ = s - (m₁ - m₂)² （入射粒子质量差的阈值）
- [1/(4Nc²) Σ|M|²] 已由 ScatteringAmplitude.jl 计算

公式参考 doc/formula/微分散射截面by散射矩阵元.md

## 设计原则

- **解耦设计**: 核心函数接受预计算的运动学变量，避免与其他模块耦合
- **高性能**: 用户可在外部预计算并复用 Mandelstam 变量和矩阵元
- **可组合性**: 适合与积分器等模块组合使用
- **运动学检查**: 提供独立的阈值和边界检查函数

## Dual Interface Pattern

This module supports **both struct and NamedTuple parameters**:

```julia
# Using structs (recommended)
using Main.ParameterTypes: QuarkParams, ThermoParams

q = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
t = ThermoParams(0.15, 0.5, 0.5, 0.0)
dσ_dt = differential_cross_section(:uu_to_uu, 2.0, -0.5, q, t, K_coeffs)

# Using NamedTuples (backward compatible)
q_nt = (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
t_nt = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
dσ_dt = differential_cross_section(:uu_to_uu, 2.0, -0.5, q_nt, t_nt, K_coeffs)
```

Both produce identical results. Internal normalization ensures type stability.

## 使用示例

```julia
# 使用 QuarkParams 和 ThermoParams 结构体（推荐）
using Main.ParameterTypes: QuarkParams, ThermoParams

q = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
t = ThermoParams(0.15, 0.5, 0.5, 0.0)

# 预计算 Mandelstam 变量
m1, m2, m3, m4 = get_quark_masses_for_process(:uu_to_uu, q)
u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
mandelstam_vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)

# 计算散射矩阵元（接受结构体参数）
M_squared = scattering_amplitude_squared(
    :uu_to_uu, s, t, q, t, K_coeffs
)

# 计算微分截面（核心函数）
dsigma_dt = differential_cross_section(
    mandelstam_vars.s_12_plus, mandelstam_vars.s_12_minus, M_squared
)

# 使用 NamedTuple（向后兼容）
q_nt = (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
t_nt = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
# ... 其余代码相同
```
"""

# Include-once helper — REMOVED (loaded via RelaxTime.jl entry point)

const _KINEMATIC_CHECKS_PATH = normpath(joinpath(@__DIR__, "KinematicChecks.jl"))
if !isdefined(Main, :KinematicChecks)
    Base.include(Main, _KINEMATIC_CHECKS_PATH)
end
using ..KinematicChecks: check_kinematic_threshold, emit_regularization_notice

# 导入参数类型
using Main.ParameterTypes: QuarkParams, ThermoParams
using Main.ParameterAdapters: normalize_quark_input, normalize_thermo_input

# 预计算常数因子：1/(16π)
const KINEMATIC_PREFACTOR = 1.0 / (16π)

# 数值容差
const EPS_REGULARIZATION = 1e-14  # s_12_minus 正则化容差

export differential_cross_section
export check_kinematic_threshold

# ----------------------------------------------------------------------------
# 核心函数：微分散射截面计算
# ----------------------------------------------------------------------------

"""
    differential_cross_section(s_12_plus, s_12_minus, M_squared) -> Float64

从预计算的运动学变量和散射矩阵元计算微分散射截面。

# 物理公式
dσ/dt = [1/(16π s₁₂⁺ s₁₂⁻)] · |M|²

其中 |M|² 已包含色-自旋平均因子 1/(4Nc²)。

# 参数
- `s_12_plus::Float64`: s - (m₁ + m₂)² [fm⁻²]
- `s_12_minus::Float64`: s - (m₁ - m₂)² [fm⁻²]  
- `M_squared::Float64`: 散射矩阵元平方 |M|² [fm⁻⁴]

# 返回值
- `Float64`: 微分散射截面 dσ/dt [fm²]

# 运动学约束
- `s_12_plus > 0`: 入射能量必须超过产生两粒子的阈值
- `s_12_minus ≠ 0`: 当 \$m_1 = m_2\$ 时需要特殊处理（自动正则化）

# 错误处理
- 如果 `s_12_plus ≤ 0`，抛出错误（违反阈值条件）
- 如果 `|s_12_minus| < 1e-14`，自动正则化并发出警告

# 参数类型支持
本函数接受预计算的运动学变量，因此不直接使用 `QuarkParams` 或 `ThermoParams`。
然而，这些参数类型在上游计算中使用（通过 `get_quark_masses_for_process` 和 
`scattering_amplitude_squared`），本模块完全支持基于结构体的工作流程。

# 示例
```julia
# 使用 QuarkParams 结构体（推荐）
using Main.ParameterTypes: QuarkParams, ThermoParams

q = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
t = ThermoParams(0.15, 0.5, 0.5, 0.0)

# 从结构体获取质量
m1, m2, m3, m4 = get_quark_masses_for_process(:uu_to_uu, q)
# 计算运动学变量
s_plus = 20.5   # fm⁻²
s_minus = 15.3  # fm⁻²
# 计算散射矩阵元（接受结构体）
M_sq = scattering_amplitude_squared(:uu_to_uu, s, t_var, q, t, K_coeffs)

dsigma_dt = differential_cross_section(s_plus, s_minus, M_sq)
println("dσ/dt = ", dsigma_dt, " fm²")

# 使用 NamedTuple（向后兼容）
q_nt = (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
# ... 其余代码相同
```

# 参考
- 公式文档: doc/formula/微分散射截面by散射矩阵元.md
- 散射矩阵元: ScatteringAmplitude.scattering_amplitude_squared
- 参数类型: Main.ParameterTypes (QuarkParams, ThermoParams)
"""
function differential_cross_section(
    s_12_plus::Float64,
    s_12_minus::Float64,
    M_squared::Float64
)::Float64
    # 检查运动学阈值
    if s_12_plus <= 0.0
        error("Kinematic threshold violation: s_12_plus = $s_12_plus ≤ 0. " *
              "This indicates s < (m1 + m2)², violating energy conservation.")
    end
    
    # 处理 m1 ≈ m2 的退化情况（s_12_minus ≈ 0）
    s_12_minus_reg = s_12_minus
    if abs(s_12_minus) < EPS_REGULARIZATION
        emit_regularization_notice(
            "s_12_minus is very small; this typically occurs when m1 ≈ m2",
            s_12_minus,
            EPS_REGULARIZATION,
        )
        s_12_minus_reg = sign(s_12_minus) * EPS_REGULARIZATION
        if s_12_minus == 0.0  # 完全相等
            s_12_minus_reg = EPS_REGULARIZATION
        end
    end
    
    # 计算运动学因子：1/(16π s_12_plus s_12_minus)
    kinematic_factor = KINEMATIC_PREFACTOR / (s_12_plus * s_12_minus_reg)
    
    # 计算微分截面
    dsigma_dt = kinematic_factor * M_squared
    
    return dsigma_dt
end

# ----------------------------------------------------------------------------
# 辅助函数：运动学约束检查
# ----------------------------------------------------------------------------

end  # module DifferentialCrossSection
