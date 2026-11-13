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

## 使用示例

本模块只接受预计算的三个数值；`QuarkParams` / `ThermoParams` 的结构体或
NamedTuple 支持属于上游 `ScatteringAmplitude.scattering_amplitude_squared`。

```julia
s = 31.0
m1 = m2 = 1.52
s_plus = s - (m1 + m2)^2
s_minus = s - (m1 - m2)^2
M_squared = 100.0  # 示例：无量纲的色自旋平均矩阵元平方

dsigma_dt = differential_cross_section(s_plus, s_minus, M_squared)  # fm⁴
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
- `M_squared::Float64`: 散射矩阵元平方 |M|² [无量纲]

# 返回值
- `Float64`: 微分散射截面 dσ/dt [fm⁴]

# 运动学约束
- `s_12_plus > 0`: 入射能量必须超过产生两粒子的阈值
- `s_12_minus ≠ 0`: 等质量时该量等于 s；接近零的输入由数值正则化处理

# 错误处理
- 如果 `s_12_plus ≤ 0`，抛出错误（违反阈值条件）
- 如果 `|s_12_minus| < 1e-14`，自动正则化并发出警告

# 参数类型支持
本函数接受预计算的运动学变量，因此不直接使用 `QuarkParams` 或 `ThermoParams`。
然而，这些参数类型在上游计算中使用（通过 `get_quark_masses_for_process` 和 
`scattering_amplitude_squared`），本模块完全支持基于结构体的工作流程。

# 示例
```julia
s = 31.0
m = 1.52
M_squared = 100.0  # 无量纲，示例预计算值
dsigma_dt = differential_cross_section(s - 4m^2, s, M_squared)
println("dσ/dt = ", dsigma_dt, " fm⁴")
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
    
    # 对接近零的运动学分母作数值正则化
    s_12_minus_reg = s_12_minus
    if abs(s_12_minus) < EPS_REGULARIZATION
        emit_regularization_notice(
            "s_12_minus is very small; applying denominator regularization",
            s_12_minus,
            EPS_REGULARIZATION,
        )
        s_12_minus_reg = sign(s_12_minus) * EPS_REGULARIZATION
        if s_12_minus == 0.0
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
