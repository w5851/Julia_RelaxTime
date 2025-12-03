# TotalPropagator API 文档

## 模块概述

`TotalPropagator` 模块用于计算PNJL模型中夸克弹性散射的总传播子。模块将多种介子（π, K, η/η', σ_π, σ_K, σ/σ'）的传播子按味因子加权求和，并支持按标量(S)和赝标量(P)通道分离计算。

**核心功能：**
- 支持所有11种散射过程（4种qq散射 + 7种qqbar散射）
- 自动识别散射道（t/u/s）对应的味因子
- 灵活处理一般介子和混合介子的传播子计算
- 支持通道分离（S/P通道），用于散射矩阵元计算
- 内部调用带缓存的极化函数，自动复用相同参数的计算结果
- 提供质心系动量计算功能

**公式参考：** `doc/formula/总传播子byPropagator.md`

---

## 主要接口函数

### 1. 通道分离接口（推荐用于散射矩阵元计算）

#### `calculate_all_propagators_by_channel`

按标量(S)和赝标量(P)通道分离计算所有散射道的总传播子。

```julia
calculate_all_propagators_by_channel(process::Symbol, k0::Float64, k_norm::Float64,
                                    quark_params::NamedTuple, thermo_params::NamedTuple,
                                    K_coeffs::NamedTuple) -> NamedTuple
```

**参数：**

| 参数 | 类型 | 单位 | 说明 |
|------|------|------|------|
| `process` | `Symbol` | - | 散射过程符号（如 `:uu_to_uu`, `:uubar_to_uubar`） |
| `k0` | `Float64` | fm⁻¹ | 介子四动量能量分量 |
| `k_norm` | `Float64` | fm⁻¹ | 介子三动量大小 |
| `quark_params` | `NamedTuple` | - | 夸克参数（见下方说明） |
| `thermo_params` | `NamedTuple` | - | 热力学参数（见下方说明） |
| `K_coeffs` | `NamedTuple` | - | 有效耦合常数系数 |

**`quark_params` 结构：**
```julia
(m = (u=m_u, d=m_d, s=m_s),   # 夸克质量（fm⁻¹）
 μ = (u=μ_u, d=μ_d, s=μ_s),   # 夸克化学势（fm⁻¹）
 A = (u=A_u, d=A_u, s=A_s))   # A函数值（fm）
```

**`thermo_params` 结构：**
```julia
(T=T, Φ=Φ, Φbar=Φbar, ξ=ξ)  # 温度（fm⁻¹）、Polyakov环、各向异性参数
```

**返回值：**

- **qq散射**：`(t_S, t_P, u_S, u_P)`，单位：fm²
  - `t_S`: t道标量通道传播子
  - `t_P`: t道赝标量通道传播子
  - `u_S`: u道标量通道传播子
  - `u_P`: u道赝标量通道传播子

- **qqbar散射**：`(t_S, t_P, s_S, s_P)`，单位：fm²
  - `t_S`: t道标量通道传播子
  - `t_P`: t道赝标量通道传播子
  - `s_S`: s道标量通道传播子
  - `s_P`: s道赝标量通道传播子

**通道分离规则：**
- **P通道（赝标量）**：π + K + η/η' 混合介子
- **S通道（标量）**：σ_π + σ_K + σ/σ' 混合介子
- 每个介子仅在其物理对应的通道中计算，无跨通道贡献

**示例：**
```julia
# qq散射
result = calculate_all_propagators_by_channel(
    :uu_to_uu, k0, k_norm, quark_params, thermo_params, K_coeffs
)
D_t_S = result.t_S  # t道标量传播子
D_t_P = result.t_P  # t道赝标量传播子

# qqbar散射
result = calculate_all_propagators_by_channel(
    :uubar_to_ssbar, k0, k_norm, quark_params, thermo_params, K_coeffs
)
D_s_S = result.s_S  # s道标量传播子
D_s_P = result.s_P  # s道赝标量传播子
```

---

### 2. 标准接口（不分离通道）

#### `calculate_all_propagators`

自动计算指定散射过程所有相关散射道的总传播子（S+P通道总和）。

```julia
calculate_all_propagators(process::Symbol, k0::Float64, k_norm::Float64,
                          quark_params::NamedTuple, thermo_params::NamedTuple,
                          K_coeffs::NamedTuple) -> NamedTuple
```

**参数：** 同 `calculate_all_propagators_by_channel`

**返回值：**

- **qq散射**：`(t = D_t, u = D_u)`，单位：fm²
- **qqbar散射**：`(t = D_t, s = D_s)`，单位：fm²

**示例：**
```julia
result = calculate_all_propagators(
    :uu_to_uu, k0, k_norm, quark_params, thermo_params, K_coeffs
)
D_t_total = result.t  # t道总传播子（S+P）
D_u_total = result.u  # u道总传播子（S+P）
```

---

### 3. 质心系动量计算

#### `calculate_cms_momentum`

根据Mandelstam变量计算质心系下介子的四动量 (k0, k)。

```julia
calculate_cms_momentum(process::Symbol, s::Float64, t::Float64, channel::Symbol,
                      quark_params::NamedTuple; u::Union{Float64, Nothing}=nothing) 
                      -> NamedTuple
```

**参数：**

| 参数 | 类型 | 单位 | 说明 |
|------|------|------|------|
| `process` | `Symbol` | - | 散射过程符号 |
| `s` | `Float64` | fm⁻² | Mandelstam变量 s |
| `t` | `Float64` | fm⁻² | Mandelstam变量 t |
| `channel` | `Symbol` | - | 散射道类型（`:s`, `:t`, `:u`） |
| `quark_params` | `NamedTuple` | - | 夸克参数（至少需要 `m` 字段） |
| `u` (可选) | `Float64` | fm⁻² | Mandelstam变量 u（默认自动计算） |

**返回值：**
```julia
(k0 = k0, k = k)
```
- `k0`: 介子能量分量（fm⁻¹）
- `k`: 介子三动量大小（fm⁻¹）

**计算规则：**
1. **s道**（q₁+q₂→介子→q₃+q₄）：k0 = |E1+E2|, k = 0
2. **t道**（q₁→q₃+介子）：k0 = |E1-E3|, k = √(k0²-t)
3. **u道**（q₁→q₄+介子）：k0 = |E1-E4|, k = √(k0²-u)

其中 E_i 为四个粒子在质心系的能量，通过 s 和质量计算。

**Mandelstam约束：**
s + t + u = m₁² + m₂² + m₃² + m₄²

如果不提供 `u` 参数，函数将自动计算以满足该约束。

**边界条件处理：**
- 当 k0²-t/u ≥ 0：正常计算 k = √(k0²-t/u)
- 当 -1e-12 < k0²-t/u < 0：数值误差，设 k = 0（发出警告）
- 当 k0²-t/u < -1e-12：明显的数值问题，设 k = 0（发出警告）

**示例：**
```julia
quark_params = (m = (u=0.3, d=0.3, s=0.5),)
s = 4.0  # fm⁻²
t = -0.5  # fm⁻²

# 自动计算 u
result = calculate_cms_momentum(:uu_to_uu, s, t, :t, quark_params)
k0 = result.k0
k = result.k

# 手动提供 u
u = 4 * 0.3^2 - s - t
result = calculate_cms_momentum(:uu_to_uu, s, t, :u, quark_params; u=u)
```

---

### 4. 质量提取辅助函数

#### `get_quark_masses_for_process`

根据散射过程符号提取四个粒子的质量。

```julia
get_quark_masses_for_process(process::Symbol, quark_params::NamedTuple) 
    -> NTuple{4, Float64}
```

**参数：**

| 参数 | 类型 | 说明 |
|------|------|------|
| `process` | `Symbol` | 散射过程符号（如 `:uu_to_uu`, `:uubar_to_ssbar`） |
| `quark_params` | `NamedTuple` | 夸克参数（需包含 `m = (u=..., d=..., s=...)` 字段） |

**返回值：**

返回元组 `(m1, m2, m3, m4)`，对应散射过程 q₁+q₂→q₃+q₄ 的四个粒子质量（单位：fm⁻¹）。

**质量映射规则：**
- 夸克和反夸克使用相同的质量（如 u 和 ū 都用 m_u）
- 从 process 符号解析粒子类型，映射到 `quark_params.m` 中的对应质量

**示例：**
```julia
quark_params = (m = (u=0.3, d=0.3, s=0.5),)
m1, m2, m3, m4 = get_quark_masses_for_process(:uubar_to_ssbar, quark_params)
# 返回: (0.3, 0.3, 0.5, 0.5)  # (m_u, m_u, m_s, m_s)
```

---

### 5. 低级接口（手动控制）

#### `total_propagator_simple`

计算一般介子（π, K, σ_π, σ_K）的总传播子，包含味因子。

```julia
total_propagator_simple(process::Symbol, channel::Symbol, meson_list::Vector{Symbol},
                       k0::Float64, k_norm::Float64,
                       quark_params::NamedTuple, thermo_params::NamedTuple,
                       K_coeffs::NamedTuple) -> ComplexF64
```

**参数：**

| 参数 | 类型 | 单位 | 说明 |
|------|------|------|------|
| `process` | `Symbol` | - | 散射过程符号 |
| `channel` | `Symbol` | - | 散射道类型（`:t`, `:u`, `:s`） |
| `meson_list` | `Vector{Symbol}` | - | 介子列表（如 `[:pi, :K]`） |
| 其他参数 | 同上 | - | - |

**返回值：** 总传播子 `ComplexF64`（单位：fm²）

**公式：**
```
D_total = T1 * sum(D_meson) * T2
```
其中 T1, T2 为味因子，通过散射过程和散射道自动确定。

---

#### `total_propagator_mixed`

计算混合介子（η/η' 或 σ/σ'）的总传播子。

```julia
total_propagator_mixed(process::Symbol, channel::Symbol, meson_channel::Symbol,
                      k0::Float64, k_norm::Float64,
                      quark_params::NamedTuple, thermo_params::NamedTuple,
                      K_coeffs::NamedTuple) -> ComplexF64
```

**参数：**

| 参数 | 类型 | 说明 |
|------|------|------|
| `meson_channel` | `Symbol` | 介子通道类型（`:P` 为 η/η'，`:S` 为 σ/σ'） |
| 其他参数 | 同上 | - |

**返回值：** 混合介子传播子 `ComplexF64`（单位：fm²）

**注意：** 混合介子传播子已包含味结构，无需额外乘味因子。

---

## 支持的散射过程

模块支持以下11种散射过程（定义在 `Constants_PNJL.SCATTERING_MESON_MAP`）：

### qq散射（4种）
| 符号 | 散射过程 | 可用散射道 |
|------|----------|-----------|
| `:uu_to_uu` | u + u → u + u | t, u |
| `:ss_to_ss` | s + s → s + s | t, u |
| `:ud_to_ud` | u + d → u + d | t, u |
| `:us_to_us` | u + s → u + s | t, u |

### qqbar散射（7种）
| 符号 | 散射过程 | 可用散射道 |
|------|----------|-----------|
| `:udbar_to_udbar` | u + đ → u + đ | t, s |
| `:usbar_to_usbar` | u + š → u + š | t, s |
| `:uubar_to_uubar` | u + ū → u + ū | t, s |
| `:uubar_to_ddbar` | u + ū → d + đ | t, s |
| `:uubar_to_ssbar` | u + ū → s + š | t, s |
| `:ssbar_to_uubar` | s + š → u + ū | t, s |
| `:ssbar_to_ssbar` | s + š → s + š | t, s |

---

## 介子类型与通道对应关系

| 介子 | 类型 | 通道 | 极化函数参数 |
|------|------|------|-------------|
| π (pi) | 一般 | P | (m_u, m_u, μ_u, μ_u, A_u, A_u) |
| K | 一般 | P | (m_u, m_s, μ_u, μ_s, A_u, A_s) |
| σ_π (sigma_pi) | 一般 | S | (m_u, m_u, μ_u, μ_u, A_u, A_u) |
| σ_K (sigma_K) | 一般 | S | (m_u, m_s, μ_u, μ_s, A_u, A_s) |
| η/η' | 混合 | P | Π_uu, Π_ss + 耦合矩阵 |
| σ/σ' | 混合 | S | Π_uu, Π_ss + 耦合矩阵 |

**说明：**
- **一般介子**：通过 `meson_propagator_simple` 计算，需要单个极化函数
- **混合介子**：通过 `meson_propagator_mixed` 计算，需要 Π_uu 和 Π_ss 两个极化函数及耦合矩阵

---

## 味因子表（表5.3）

散射道的味因子 T1 和 T2 通过以下规则确定：

- **t道**：q1→q3 (T1), q2→q4 (T2)
- **u道**：q1→q4 (T1), q2→q3 (T2)
- **s道**：q1-q2 (T1), q3-q4 (T2)

基础味因子：

| (flavor1, flavor2) | 味因子值 |
|-------------------|---------|
| (u, u) | 1.0 |
| (d, d) | -1.0 |
| (s, s) | 2.0 |
| (u, d) | √2 |
| (u, s) | √2 |
| (d, s) | √2 |

**对称性：** T(q1, q2) = T(q2, q1)

---

## 缓存管理

模块内部使用 `PolarizationCache` 自动缓存极化函数的计算结果。

### 缓存控制函数

```julia
reset_cache!()  # 清空缓存
get_cache_stats()  # 获取缓存统计信息
```

**使用建议：**
- 在温度扫描等批量计算前调用 `reset_cache!()` 清空缓存
- 使用 `get_cache_stats()` 监控缓存命中率

---

## 完整示例

### 示例1：使用通道分离接口计算qq散射

```julia
using .TotalPropagator
using .EffectiveCouplings
using .OneLoopIntegrals: A
using .Constants_PNJL

# 1. 设置物理参数
T = 150.0 / 197.327  # 150 MeV → fm⁻¹
m_u = 300.0 / 197.327
m_s = 500.0 / 197.327
μ_u = 0.0
μ_s = 0.0
Φ = 0.5
Φbar = 0.5
ξ = 0.0  # 各向同性

# 2. 计算夸克参数
A_u = A(T, μ_u, m_u, Φ, Φbar)
A_s = A(T, μ_s, m_s, Φ, Φbar)

G_u = calculate_G_from_A(A_u)
G_s = calculate_G_from_A(A_s)

quark_params = (
    m = (u=m_u, d=m_u, s=m_s),
    μ = (u=μ_u, d=μ_u, s=μ_s),
    A = (u=A_u, d=A_u, s=A_s)
)

thermo_params = (T=T, Φ=Φ, Φbar=Φbar, ξ=ξ)

# 3. 计算有效耦合常数
G_fm2 = Constants_PNJL.G_GeV_inv2 / (197.327^2)
K_fm5 = Constants_PNJL.K_GeV_inv5 / (197.327^5)
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

# 4. 计算质心系动量
s = 4.0  # fm⁻²
t = -0.5  # fm⁻²
result_t = calculate_cms_momentum(:uu_to_uu, s, t, :t, quark_params)
k0_t = result_t.k0
k_t = result_t.k

result_u = calculate_cms_momentum(:uu_to_uu, s, t, :u, quark_params)
k0_u = result_u.k0
k_u = result_u.k

# 5. 计算通道分离的传播子
# 方法A：分别计算t道和u道
D_t_result = calculate_all_propagators_by_channel(
    :uu_to_uu, k0_t, k_t, quark_params, thermo_params, K_coeffs
)
println("t道标量传播子: ", D_t_result.t_S)
println("t道赝标量传播子: ", D_t_result.t_P)

# 方法B：如果k0和k对两个道相同（简化情况）
# result = calculate_all_propagators_by_channel(
#     :uu_to_uu, k0, k_norm, quark_params, thermo_params, K_coeffs
# )
# println("t道: S=", result.t_S, ", P=", result.t_P)
# println("u道: S=", result.u_S, ", P=", result.u_P)
```

### 示例2：使用标准接口（不分离通道）

```julia
# 使用相同的参数设置

k0 = 100.0 / 197.327
k_norm = 50.0 / 197.327

# 计算总传播子（S+P）
result = calculate_all_propagators(
    :uu_to_uu, k0, k_norm, quark_params, thermo_params, K_coeffs
)

println("t道总传播子: ", result.t)
println("u道总传播子: ", result.u)

# 验证一致性（分离 vs 不分离）
result_sep = calculate_all_propagators_by_channel(
    :uu_to_uu, k0, k_norm, quark_params, thermo_params, K_coeffs
)
@assert result.t ≈ result_sep.t_S + result_sep.t_P
@assert result.u ≈ result_sep.u_S + result_sep.u_P
```

### 示例3：qqbar散射

```julia
# 计算 u + ū → s + š 过程的传播子
s = 5.0  # fm⁻²
t = -0.8  # fm⁻²

# 计算s道和t道的质心系动量
result_s = calculate_cms_momentum(:uubar_to_ssbar, s, t, :s, quark_params)
result_t = calculate_cms_momentum(:uubar_to_ssbar, s, t, :t, quark_params)

# 使用s道动量计算s道传播子
result_s_prop = calculate_all_propagators_by_channel(
    :uubar_to_ssbar, result_s.k0, result_s.k, 
    quark_params, thermo_params, K_coeffs
)
println("s道标量传播子: ", result_s_prop.s_S)
println("s道赝标量传播子: ", result_s_prop.s_P)

# 使用t道动量计算t道传播子
result_t_prop = calculate_all_propagators_by_channel(
    :uubar_to_ssbar, result_t.k0, result_t.k, 
    quark_params, thermo_params, K_coeffs
)
println("t道标量传播子: ", result_t_prop.t_S)
println("t道赝标量传播子: ", result_t_prop.t_P)
```

---

## 错误处理

### 常见错误及解决方法

1. **未知散射过程**
   ```julia
   ERROR: Unknown scattering process: :invalid_process. 
          Supported processes: :uu_to_uu, :ss_to_ss, ...
   ```
   **解决：** 检查散射过程符号，参考"支持的散射过程"章节。

2. **无效散射道**
   ```julia
   ERROR: Scattering process :uu_to_uu (qq type) only has t and u channels, not :s
   ```
   **解决：** qq散射只有t和u道，qqbar散射只有t和s道。

3. **k0²-t/u < 0 警告**
   ```julia
   WARNING: calculate_cms_momentum: k0²-t = -0.001 < -1e-12, setting k=0. 
            Check Mandelstam variables.
   ```
   **原因：** Mandelstam变量不满足物理约束或数值计算误差。
   **解决：** 检查输入的 s, t, u 是否满足约束 s+t+u=Σm²，或增加数值精度。

4. **无效介子类型**
   ```julia
   ERROR: Unknown meson type: :invalid_meson. Use :pi, :K, :sigma_pi, or :sigma_K
   ```
   **解决：** 使用正确的介子符号。

---

## 性能建议

1. **利用缓存**：相同参数的极化函数只计算一次，自动复用缓存结果。
2. **批量计算**：在温度扫描等场景，相同 (k0, k) 的传播子会被缓存。
3. **避免频繁清空缓存**：只在参数体系完全改变时调用 `reset_cache!()`。
4. **通道分离与否的选择**：
   - 需要S/P分量 → 使用 `calculate_all_propagators_by_channel`
   - 只需总传播子 → 使用 `calculate_all_propagators`（略快）

---

## 依赖模块

- `Constants_PNJL`：散射过程与介子映射表
- `MesonPropagator`：单个介子传播子计算
- `PolarizationCache`：极化函数缓存系统
- `EffectiveCouplings`：有效耦合常数计算
- `OneLoopIntegrals`：单圈积分（A函数）

---

## 版本历史

- **v2.0**（当前）：
  - 新增 `calculate_all_propagators_by_channel` 函数，支持S/P通道分离
  - 新增 `calculate_cms_momentum` 函数，从Mandelstam变量计算介子四动量
  - 新增 `get_quark_masses_for_process` 辅助函数
  - 扩展支持所有11种散射过程
  - 改进边界条件处理（k0²-t/u < 0的情况）

- **v1.0**：
  - 初始版本，支持 `calculate_all_propagators` 和低级接口
  - 支持一般介子和混合介子的总传播子计算
  - 集成极化函数缓存系统

---

## 参考文献

1. 公式文档：`doc/formula/总传播子byPropagator.md`
2. 散射矩阵元：`doc/formula/ScatteringAmplitude_散射矩阵元by总传播子.md`
3. 介子传播子：`api/MesonPropagator.md`
4. 极化函数：`api/PolarizationAniso.md`
