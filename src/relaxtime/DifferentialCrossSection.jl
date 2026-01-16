module DifferentialCrossSection

"""
# DifferentialCrossSection.jl

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

```julia
# 预计算 Mandelstam 变量
m1, m2, m3, m4 = get_quark_masses_for_process(:uu_to_uu, quark_params)
u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
mandelstam_vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)

# 计算散射矩阵元
M_squared = scattering_amplitude_squared(
    :uu_to_uu, s, t, quark_params, thermo_params, K_coeffs
)

# 计算微分截面（核心函数）
dsigma_dt = differential_cross_section(
    mandelstam_vars.s_12_plus, mandelstam_vars.s_12_minus, M_squared
)
```
"""

# 预计算常数因子：1/(16π)
const KINEMATIC_PREFACTOR = 1.0 / (16π)

# 数值容差
const EPS_THRESHOLD = 1e-12  # 阈值检查容差
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

# 示例
```julia
# 从 calculate_mandelstam_variables 获取
s_plus = 20.5   # fm⁻²
s_minus = 15.3  # fm⁻²
M_sq = 2468.5   # fm⁻⁴

dsigma_dt = differential_cross_section(s_plus, s_minus, M_sq)
println("dσ/dt = ", dsigma_dt, " fm²")
```

# 参考
- 公式文档: doc/formula/微分散射截面by散射矩阵元.md
- 散射矩阵元: ScatteringAmplitude.scattering_amplitude_squared
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
        @warn "s_12_minus is very small (|s_12_minus| = $(abs(s_12_minus)) < $EPS_REGULARIZATION), " *
              "applying regularization. This typically occurs when m1 ≈ m2."
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

"""
    check_kinematic_threshold(s, m1, m2; warn_close=true) -> Bool

检查 Mandelstam 变量 s 是否满足运动学阈值条件。

# 阈值条件
对于散射过程 q₁ + q₂ → q₃ + q₄，要求：
s ≥ (m₁ + m₂)²

这确保了质心系中入射粒子的相对动量为实数。

# 参数
- `s::Float64`: Mandelstam 变量 s [fm⁻²]
- `m1::Float64`: 第一个入射粒子质量 [fm⁻¹]
- `m2::Float64`: 第二个入射粒子质量 [fm⁻¹]
- `warn_close::Bool=true`: 是否在接近阈值时发出警告

# 返回值
- `true`: 满足阈值条件
- `false`: 违反阈值条件（s < (m1+m2)²）

# 警告
当 `s` 非常接近阈值时（`s_12_plus < 1e-12`），会发出警告，
因为此时微分截面可能发散。

# 示例
```julia
s = 4.0  # fm⁻²
m_u = 0.3  # fm⁻¹

if check_kinematic_threshold(s, m_u, m_u)
    println("运动学条件满足，可以计算散射截面")
else
    println("警告：s 低于阈值！")
end
```

# 参考
- 阈值物理: s_threshold = (m1 + m2)² 对应质心系零动量
- 接近阈值时的发散行为与粒子产生阈值相关
"""
function check_kinematic_threshold(
    s::Float64,
    m1::Float64,
    m2::Float64;
    warn_close::Bool=true
)::Bool
    s_threshold = (m1 + m2)^2
    
    # 检查是否低于阈值
    if s < s_threshold
        @warn "Kinematic threshold violation" s=s threshold=s_threshold deficit=(s_threshold - s)
        return false
    end
    
    # 检查是否非常接近阈值
    s_plus = s - s_threshold
    if warn_close && s_plus < EPS_THRESHOLD
        @warn "s is very close to threshold (s_12_plus = $s_plus < $EPS_THRESHOLD). " *
              "Differential cross section may diverge near threshold."
    end
    
    return true
end


end  # module DifferentialCrossSection
