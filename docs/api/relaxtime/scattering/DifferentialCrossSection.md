# DifferentialCrossSection 模块 API 文档

## 模块概述

`DifferentialCrossSection` 模块计算相对论性 Boltzmann 动力学理论中的微分散射截面，连接散射矩阵元模块与总散射截面/弛豫时间计算。

**文件位置**: `src/relaxtime/DifferentialCrossSection.jl`

## 物理背景

### 微分散射截面公式

在相对论性 Boltzmann 动力学理论中，微分散射截面为：

$$\frac{d\sigma}{dt} = \frac{1}{16\pi s_{12}^+ s_{12}^-} \cdot \frac{1}{4N_C^2} \sum_{\text{spin, color}} |\mathcal{M}|^2$$

其中：
- $s_{12}^+ = s - (m_1 + m_2)^2$：入射粒子质量和的阈值
- $s_{12}^- = s - (m_1 - m_2)^2$：入射粒子质量差的阈值  
- $\frac{1}{4N_C^2} \sum |\mathcal{M}|^2$：已由 `ScatteringAmplitude.jl` 计算

### 运动学约束

- **阈值条件**: $s \geq (m_1 + m_2)^2$ 确保质心系动量为实数
- **t 积分边界**: $t \in [t_{\min}, 0]$，其中 $t_{\min} = -\frac{[s - (m_1+m_2)^2][s - (m_3+m_4)^2]}{4s}$

## 设计原则

1. **解耦设计**: 核心函数接受预计算的运动学变量，避免与其他模块耦合
2. **高性能**: 用户可在外部预计算并复用 Mandelstam 变量和矩阵元
3. **可组合性**: 适合与积分器等模块组合使用
4. **运动学检查**: 提供独立的阈值和边界检查函数

---

## 核心函数

### `differential_cross_section`

```julia
differential_cross_section(s_12_plus, s_12_minus, M_squared) -> Float64
```

从预计算的运动学变量和散射矩阵元计算微分散射截面。

#### 参数

- `s_12_plus::Float64`: $s - (m_1 + m_2)^2$ [fm⁻²]
- `s_12_minus::Float64`: $s - (m_1 - m_2)^2$ [fm⁻²]
- `M_squared::Float64`: 散射矩阵元平方 $|\mathcal{M}|^2$ [fm⁻⁴]

#### 返回值

- `Float64`: 微分散射截面 $d\sigma/dt$ [fm²]

#### 运动学约束

- `s_12_plus > 0`: 必须超过阈值，否则抛出错误
- `|s_12_minus| > 0`: 当 $m_1 = m_2$ 时自动正则化（$\epsilon = 10^{-14}$）

#### 示例

```julia
using DifferentialCrossSection
using ScatteringAmplitude

# 预计算 Mandelstam 变量
s = 31.0  # fm⁻²
t = -2.0
m1, m2, m3, m4 = get_quark_masses_for_process(:uu_to_uu, quark_params)
u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
mandelstam_vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)

# 计算散射矩阵元
M_squared = scattering_amplitude_squared(
    :uu_to_uu, s, t, quark_params, thermo_params, K_coeffs
)

# 计算微分截面
dsigma_dt = differential_cross_section(
    mandelstam_vars.s_12_plus,
    mandelstam_vars.s_12_minus,
    M_squared
)

println("dσ/dt = ", dsigma_dt, " fm²")
# 输出: dσ/dt = 0.065074 fm²
```

#### 错误处理

- **阈值违反**: 如果 `s_12_plus ≤ 0`，抛出 `ErrorException`
- **退化情况**: 如果 `|s_12_minus| < 1e-14`，自动正则化并发出 `@warn`

---

## 辅助函数

### `check_kinematic_threshold`

```julia
check_kinematic_threshold(s, m1, m2; warn_close=true) -> Bool
```

检查 Mandelstam 变量 $s$ 是否满足运动学阈值条件。

#### 参数

- `s::Float64`: Mandelstam 变量 s [fm⁻²]
- `m1::Float64`: 第一个入射粒子质量 [fm⁻¹]
- `m2::Float64`: 第二个入射粒子质量 [fm⁻¹]
- `warn_close::Bool=true`: 是否在接近阈值时发出警告

#### 返回值

- `true`: 满足阈值条件 ($s \geq (m_1 + m_2)^2$)
- `false`: 违反阈值条件

#### 示例

```julia
s = 10.0  # fm⁻²
m_u = 1.52  # fm⁻¹

if check_kinematic_threshold(s, m_u, m_u)
    println("运动学条件满足")
else
    error("s 低于阈值，无法计算散射截面")
end

# 接近阈值时的警告
s_threshold = (m_u + m_u)^2
s_close = s_threshold + 1e-13
check_kinematic_threshold(s_close, m_u, m_u)
# 警告: s is very close to threshold...
```

---

### 注意：`calculate_t_bounds` 已移至 `TotalCrossSection` 模块

`calculate_t_bounds` 函数已从 `DifferentialCrossSection` 模块移除，现在位于 `TotalCrossSection` 模块中。

**原因**：
- t 积分边界与总散射截面计算直接相关
- 微分散射截面模块只处理 dσ/dt 的计算，不涉及积分逻辑
- 更好的模块化和职责分离

**使用方式**：

```julia
using .TotalCrossSection: calculate_t_bounds

t_bounds = calculate_t_bounds(s, m1, m2, m3, m4)
# 返回: (t_min = ..., t_max = ...)
```

**详细文档**：请查看 `api/TotalCrossSection.md`

---

## 使用场景

### 场景1：单点计算

适用于计算特定动力学点的微分截面。

```julia
# 1. 设置参数
s = 31.0
t = -2.0
process = :uu_to_uu

# 2. 预计算运动学变量
m1, m2, m3, m4 = get_quark_masses_for_process(process, quark_params)
u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
mandelstam_vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)

# 3. 计算矩阵元
M_squared = scattering_amplitude_squared(
    process, s, t, quark_params, thermo_params, K_coeffs
)

# 4. 计算微分截面
dsigma_dt = differential_cross_section(
    mandelstam_vars.s_12_plus,
    mandelstam_vars.s_12_minus,
    M_squared
)
```

### 场景2：t 积分（总截面）

适用于计算总散射截面 $\sigma_{\text{total}} = \int_{t_{\min}}^{t_{\max}} \frac{d\sigma}{dt} dt$。

**注意**：实际上应使用 `TotalCrossSection` 模块的 `total_cross_section` 函数，该函数已包含完整实现（含末态统计因子）。

```julia
using QuadGK
using .TotalCrossSection: calculate_t_bounds  # 从 TotalCrossSection 导入

function total_cross_section(s, process, quark_params, thermo_params, K_coeffs)
    # 获取质量和 t 边界
    m1, m2, m3, m4 = get_quark_masses_for_process(process, quark_params)
    t_bounds = calculate_t_bounds(s, m1, m2, m3, m4)
    
    # 预计算 s 相关的 Mandelstam 变量（t 独立部分）
    s_12_plus = s - (m1 + m2)^2
    s_12_minus = s - (m1 - m2)^2
    
    # 对 t 积分
    σ_total, err = quadgk(t_bounds.t_min, t_bounds.t_max) do t
        # 计算该 t 点的矩阵元
        M_squared = scattering_amplitude_squared(
            process, s, t, quark_params, thermo_params, K_coeffs
        )
        
        # 计算微分截面
        differential_cross_section(s_12_plus, s_12_minus, M_squared)
    end
    
    return σ_total
end

# 使用
σ_uu = total_cross_section(31.0, :uu_to_uu, quark_params, thermo_params, K_coeffs)
println("σ_total(uu→uu) = ", σ_uu, " fm²")
```

### 场景3：批量计算多个过程

适用于对所有散射过程进行统一计算。

```julia
function calculate_all_dsigma_dt(s, t, quark_params, thermo_params, K_coeffs)
    # 先批量计算所有矩阵元
    M_squared_all = calculate_all_scattering_amplitudes_squared(
        s, t, quark_params, thermo_params, K_coeffs
    )
    
    results = Dict{Symbol, Float64}()
    
    for (process, M_squared) in pairs(M_squared_all)
        # 获取该过程的质量
        m1, m2, m3, m4 = get_quark_masses_for_process(process, quark_params)
        
        # 计算运动学因子
        s_12_plus = s - (m1 + m2)^2
        s_12_minus = s - (m1 - m2)^2
        
        # 计算微分截面
        results[process] = differential_cross_section(
            s_12_plus, s_12_minus, M_squared
        )
    end
    
    return NamedTuple(results)
end

# 使用
all_dsigma = calculate_all_dsigma_dt(31.0, -2.0, quark_params, thermo_params, K_coeffs)
for (process, dsigma) in pairs(all_dsigma)
    println("$process: dσ/dt = $dsigma fm²")
end
```

---

## 性能考虑

### 优化策略

1. **预计算复用**: 
   - Mandelstam 变量在 t 积分循环外预计算
   - 质量在批量计算中提前获取

2. **常数预计算**:
   - 模块内置 `KINEMATIC_PREFACTOR = 1/(16π)`
   - 避免重复计算三角函数

3. **批量计算**:
   - 先调用 `calculate_all_scattering_amplitudes_squared`
   - 传播子只计算一次（最耗时部分）

### 典型耗时

| 操作 | 耗时 | 说明 |
|------|------|------|
| `differential_cross_section` | ~1 ns | 纯算术运算 |
| `calculate_mandelstam_variables` | ~10 ns | 18个变量计算 |
| `scattering_amplitude_squared` | ~50 ms | 含传播子计算 |
| 总截面 t 积分（100点） | ~5 s | 主要是矩阵元计算 |

**结论**: 微分截面计算本身开销极小，瓶颈在散射矩阵元的传播子计算。

---

## 物理验证

### 测试结果示例

```julia
# 参数: T=150 MeV, μ=0, m_u≈1.52 fm⁻¹, m_s≈2.53 fm⁻¹
# s=31.0 fm⁻², t=-2.0 fm⁻²

过程                 dσ/dt (fm²)
----------------------------------
uu→uu                0.065074
ss→ss                0.008253
udbar→udbar          0.262862
uubar→uubar          0.094940
```

### 物理一致性

- ✅ **同位旋对称**: $d\sigma/dt(uu) = d\sigma/dt(dd)$ 当 $m_u = m_d$
- ✅ **电荷共轭**: $d\sigma/dt(\text{dubar}) = d\sigma/dt(\text{udbar})$ （相对误差 < $10^{-12}$）
- ✅ **阈值行为**: 接近 $s \to (m_1+m_2)^2$ 时截面发散
- ✅ **正定性**: $d\sigma/dt > 0$ 对所有物理参数

---

## 依赖关系

### 导入模块

```julia
using .ScatteringAmplitude: scattering_amplitude_squared,
                            get_quark_masses_for_process,
                            calculate_mandelstam_variables
```

### 与其他模块的关系

```
ScatteringAmplitude (|M|²)
          ↓
DifferentialCrossSection (dσ/dt)
          ↓
TotalCrossSection (σ_total)
          ↓
RelaxationTime (τ)
```

---

## 注意事项

### 常见陷阱

1. **阈值检查**: 必须确保 $s \geq (m_1 + m_2)^2$，否则 `s_12_plus ≤ 0` 导致错误
2. **退化情况**: 当 $m_1 = m_2$ 时，`s_12_minus = s`，不影响计算
3. **单位一致性**: 输入 [fm⁻²], [fm⁻⁴] → 输出 [fm²]
4. **t 边界**: 正向散射 $t_{\max} = 0$，后向散射 $t_{\min} < 0$

### 调试建议

```julia
# 打开运动学检查
check_kinematic_threshold(s, m1, m2, warn_close=true)

# 检查 t 是否在物理范围内
using .TotalCrossSection: calculate_t_bounds  # 从 TotalCrossSection 导入
t_bounds = calculate_t_bounds(s, m1, m2, m3, m4)
@assert t_bounds.t_min ≤ t ≤ t_bounds.t_max "t out of physical range"

# 验证公式关系
kinematic_factor = 1.0 / (16π * s_12_plus * s_12_minus)
@assert dsigma_dt ≈ kinematic_factor * M_squared
```

---

## 参考文献

- **公式推导**: `docs/reference/formula/relaxtime/scattering/DifferentialCrossSection_FromScatteringAmplitude.md`
- **散射矩阵元**: `docs/api/relaxtime/scattering/ScatteringAmplitude.md`
- **测试用例**: `test_unit/test_differential_cross_section.jl`
