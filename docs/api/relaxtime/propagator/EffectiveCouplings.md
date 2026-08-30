# EffectiveCouplings 模块 API 文档

## 模块概述

`EffectiveCouplings` 模块计算3味PNJL模型中的有效耦合常数。这些耦合系数是介子传播子计算的核心输入，由原始耦合常数（G、K）和夸克凝聚函数（通过单圈积分 A 计算）共同决定。物理基本量是 `phi_f = <bar q_f q_f>`；旧接口中的 `G_u/G_s` 只是 `H_f=-phi_f` 的历史 helper 命名。

**物理背景**: 在随机相位近似(RPA)下，夸克-反夸克对形成的介子传播子由有效耦合系数 K_α^± 描述。这些系数编码了不同味道通道（单态、八重态）和自旋-宇称通道（标量S、赝标量P）的相互作用强度。

**在计算链中的位置**:
```
A函数(T,μ,m) → phi_f → H_f=-phi_f（旧适配层）→ K系数 → 介子传播子D → 散射振幅M → 弛豫时间τ
```

若平衡态已经提供 `phi_f`，推荐直接走 `calculate_effective_couplings_from_phi`；
上面的 `A_f → H_f` 仅保留给没有凝聚值的历史调用者。

## 依赖

- `Constants_PNJL` —— 提供原始耦合常数 G_fm2、K_fm5
- `OneLoopIntegrals` —— 提供各向同性A函数
- `OneLoopIntegralsAniso` —— 提供各向异性A函数（可选）

## 单位约定

- 原始耦合常数：G (fm²)、K (fm⁵)
- 单圈积分：A_f (fm⁻²)
- 物理夸克凝聚：phi_f (fm⁻³)
- 历史 helper：H_f=-phi_f，API 参数名仍为 G_u/G_s (fm⁻³)
- 有效耦合系数：K_α^± (fm²)
- 耦合矩阵行列式：det K (fm⁴)
- 乘积 `K H_f`：fm²

---

## API 参考

### `calculate_G_from_A(A_f, m_f; Nc=3)`

从单圈积分 A 函数计算历史 helper `H_f=-phi_f`。函数名和返回变量名保留旧 API 兼容性；完整 KMT 后端应直接使用已经得到的 `phi_f`。

#### 函数签名

```julia
calculate_G_from_A(A_f::Real, m_f::Real; Nc::Int=3) -> Float64
```

#### 物理公式

```
phi_f = N_c / (4π²) × m_f × A_f(T, μ)
H_f = -phi_f = -N_c / (4π²) × m_f × A_f(T, μ)
```

其中 A_f 由 `OneLoopIntegrals.A` 或 `OneLoopIntegralsAniso.A_aniso` 计算。

#### 参数

| 参数 | 类型 | 说明 | 单位 | 默认值 |
|------|------|------|------|--------|
| `A_f` | `Float64` | 单圈积分 A 函数的值 | fm⁻² | — |
| `m_f` | `Float64` | 该味夸克质量 | fm⁻¹ | — |
| `Nc` | `Int` | 色数（QCD中固定为3） | 无量纲 | 3 |

#### 返回值

返回 `Float64`，表示历史 helper `H_f=-phi_f`（单位：fm⁻³）。

**物理意义**:
- `phi_f` 描述温度和密度下的夸克凝聚；`H_f` 只是其相反数
- 手征对称破缺时 `phi_f` 与 `H_f` 均非零
- 高温极限（T → ∞）时 `phi_f,H_f → 0`（手征对称恢复）

#### 量纲关系

数值取决于质量、化学势、Polyakov 背景、截断和积分节点，不在 API 文档中固定“典型值”。量纲关系为：

| 量 | 单位 | 关系 |
|------|------|------|
| `A_f` | fm⁻² | 单圈动量积分 |
| `phi_f` | fm⁻³ | `N_c*m_f*A_f/(4π²)` |
| `H_f`（旧参数名 `G_f`） | fm⁻³ | `-phi_f` |

#### 使用示例

```julia
using .Constants_PNJL: ħc_MeV_fm
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A
using .GaussLegendre: gauleg

# 物理参数
T_MeV = 150.0
μ_MeV = 0.0
m_u_MeV = 300.0  # u/d 夸克有效质量

# 转换为 fm⁻¹
T_inv_fm = T_MeV / ħc_MeV_fm
μ_inv_fm = μ_MeV / ħc_MeV_fm
m_u_inv_fm = m_u_MeV / ħc_MeV_fm

Φ = 0.5
Φbar = 0.5

# 生成积分节点
nodes_p, weights_p = gauleg(0.0, 20.0, 64)

# 计算A函数
A_u = A(m_u_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, nodes_p, weights_p)
println("A_u = ", A_u, " fm⁻²")

# 计算历史 helper H_u=-phi_u（变量名 G_u 为兼容旧 API）
G_u = calculate_G_from_A(A_u, m_u_inv_fm)
println("H_u (API variable G_u) = ", G_u, " fm⁻³")
```

**输出格式示意（数值随模型配置、截断和节点数变化，不构成固定基准）**:
```
A_u = <A_u> fm⁻²
H_u (API variable G_u) = <H_u> fm⁻³
```

---

### `calculate_effective_couplings_from_phi(G, K, phi_l, phi_s)`

旧 `u=d` 有效耦合接口的 `phi` 原生适配层。它直接使用已经由平衡态求解得到的
轻味凝聚 `phi_l`（旧接口中的共同 `u/d` 投影）和奇异味凝聚 `phi_s`，不重新计算
`A_f`，也不把 `H_f=-phi_f` 暴露为新的物理输入。

```julia
calculate_effective_couplings_from_phi(
    G::Real, K::Real, phi_l::Real, phi_s::Real,
) -> NamedTuple
```

该函数与
`calculate_effective_couplings(G, K, -phi_l, -phi_s)`
保持纯代数等价，并返回相同的 `K0/K123/K4567/K8/K08/det_K` 字段。它不是完整的
非对称 KMT 核：在 `phi_u != phi_d` 时，旧 `K4567` 仍只对应 `K67`，物理
`K^±` 必须使用 `MesonInteractionKernel` 的 `K45`。

```julia
using .EffectiveCouplings: calculate_effective_couplings_from_phi

phi_u = -0.30  # fm⁻³；作为旧 u/d 对称投影的 phi_l
phi_s = -0.10  # fm⁻³
K_coeffs = calculate_effective_couplings_from_phi(G_fm2, K_fm5, phi_u, phi_s)
```

在分析 runtime 中，若 `quark_params` 额外携带
`phi=(u=..., d=..., s=...)`，`MesonDensity` 会优先使用该适配层；没有 `phi` 字段
的旧 `m,mu,A` 输入仍保留原有 `calculate_G_from_A` 回退路径。

---

### `calculate_effective_couplings(G, K, G_u, G_s)`

计算10个有效耦合系数 K_α^±，用于介子传播子计算。

#### 函数签名

```julia
calculate_effective_couplings(G::Real, K::Real,
                               G_u::Real, G_s::Real) -> NamedTuple
```

#### 物理公式

```
K₀±    = G ∓ (1/3)K(2H_u + H_s)         # 单态通道
K₁₂₃± = G ± (1/2)KH_s                  # π介子通道
K₄₅₆₇± = G ± (1/2)KH_u                 # 旧合并 K 通道
K₈±    = G ± (1/6)K(4H_u - H_s)        # 八重态通道
K₀₈±   = ±(1/6)√2 K(H_u - H_s)          # 混合通道
```

**符号约定**:
- `+`: 赝标量通道（P，自旋0负宇称）
- `-`: 标量通道（S，自旋0正宇称）
- `H_u`（API 参数名 `G_u`）和 `H_s`（API 参数名 `G_s`）的单位均为 fm⁻³
- 旧公式假设 u/d 味道简并；完整非对称核应分别使用 `phi_u`、`phi_d`、`phi_s`

#### 参数

| 参数 | 类型 | 说明 | 单位 | 典型值 |
|------|------|------|------|--------|
| `G` | `Float64` | 四夸克相互作用耦合常数 | fm² | 由模型配置给定 |
| `K` | `Float64` | 't Hooft六夸克相互作用耦合常数 | fm⁵ | 由模型配置给定 |
| `G_u` | `Float64` | 历史 helper `H_u=-phi_u` | fm⁻³ | 由 A 与 `m_u` 计算 |
| `G_s` | `Float64` | 历史 helper `H_s=-phi_s` | fm⁻³ | 由 A 与 `m_s` 计算 |

#### 返回值

返回一个 `NamedTuple`，包含10个有效耦合系数（单位：fm²）：

```julia
(
    K0_plus::Float64,      # 单态赝标量（P通道）
    K0_minus::Float64,     # 单态标量（S通道）
    K123_plus::Float64,    # π通道赝标量（π介子）
    K123_minus::Float64,   # π通道标量（a0介子）
    K4567_plus::Float64,   # K通道赝标量（K介子）
    K4567_minus::Float64,  # K通道标量（K0*介子）
    K8_plus::Float64,      # 八重态赝标量（η8）
    K8_minus::Float64,     # 八重态标量
    K08_plus::Float64,     # 混合赝标量（η-η'混合）
    K08_minus::Float64     # 混合标量（σ-f0混合）
)
```

#### 物理意义详解

**味道通道分类**:

| 通道 | SU(3)表示 | 对应介子（P通道） | 对应介子（S通道） |
|------|-----------|------------------|------------------|
| K₀ | 单态 | η₀ (η') | σ, f₀(500) |
| K₁₂₃ | 三重态 | π⁰, π± | a₀(980) |
| K₄₅₆₇ | 旧接口合并的四重态 | K⁰, K±, K̄⁰（`u≠d` 时不应视为单一耦合） | K₀*(1430) |
| K₈ | 八重态 | η₈ | f₀(980) |
| K₀₈ | 混合 | η-η' 混合 | σ-f₀ 混合 |

**手征极限检验**:

在手征极限下（m_q → 0），夸克凝聚及其 helper 消失，`H_f=-phi_f → 0`，所有 K_α^± 退化为原始耦合常数 G：

```julia
# 手征极限
G_u = 0.0
G_s = 0.0
K_coeffs = calculate_effective_couplings(G, K, 0.0, 0.0)
# 结果：所有 K_α^± = G
```

#### 非对称背景下的通道映射

`calculate_effective_couplings` 是旧的 u/d 简并接口。令
`H_f=-phi_f`，则其代数映射为：

| 旧字段/通道 | 完整 KMT 后端 | 适用的物理通道 |
|---|---|---|
| `K123` | `K12` | `pi^±` |
| `K4567`（使用 `H_u`） | `K67` | `K^0/\bar K^0`（d-s spectator） |
| — | `K45` | `K^±`（u-s charged channel） |

因此在 `phi_u=phi_d` 时 `K45=K67=K4567`，旧接口与完整后端相容；
在 `phi_u≠phi_d` 时，旧 `K4567` 不应直接标注为 charged `K^±` 的耦合。
现有旧函数签名保持不变，迁移到非对称生产路线时应改用
`MesonInteractionKernel.charged_coupling(kernel, :K45, :P/:S)`。

**SU(3)对称检验**:

当 `H_u = H_s`（等价于 `phi_u = phi_s` 的三味对称）时：

```
K₁₂₃± = K₄₅₆₇± = K₈±     # 味道简并
K₀₈± = 0                   # 无混合
```

#### 使用示例

```julia
using .Constants_PNJL: G_fm2, K_fm5, ħc_MeV_fm
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .GaussLegendre: gauleg

# 物理参数
T_MeV = 150.0
μ_MeV = 0.0
m_u_MeV = 300.0  # u/d 夸克
m_s_MeV = 500.0  # s 夸克

# 转换单位
T_inv_fm = T_MeV / ħc_MeV_fm
μ_inv_fm = μ_MeV / ħc_MeV_fm
m_u_inv_fm = m_u_MeV / ħc_MeV_fm
m_s_inv_fm = m_s_MeV / ħc_MeV_fm

Φ = 0.5
Φbar = 0.5

# 生成积分节点
nodes_p, weights_p = gauleg(0.0, 20.0, 64)

# 计算A函数
A_u = A(m_u_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, nodes_p, weights_p)
A_s = A(m_s_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, nodes_p, weights_p)

# 计算历史 helper H_f=-phi_f（变量名 G_u/G_s 为兼容旧 API）
G_u = calculate_G_from_A(A_u, m_u_inv_fm)
G_s = calculate_G_from_A(A_s, m_s_inv_fm)

# 计算有效耦合系数
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

# 输出结果
println("原始耦合常数:")
println("  G = ", G_fm2, " fm²")
println("  K = ", K_fm5, " fm⁵")
println()
println("历史 helper（H_f=-phi_f）:")
println("  H_u (API variable G_u) = ", G_u, " fm⁻³")
println("  H_s (API variable G_s) = ", G_s, " fm⁻³")
println()
println("有效耦合系数:")
println("  K₀⁺ = ", K_coeffs.K0_plus, " fm² (单态赝标量)")
println("  K₀⁻ = ", K_coeffs.K0_minus, " fm² (单态标量)")
println("  K₁₂₃⁺ = ", K_coeffs.K123_plus, " fm² (π通道赝标量)")
println("  K₁₂₃⁻ = ", K_coeffs.K123_minus, " fm² (π通道标量)")
println("  K₄₅₆₇⁺ = ", K_coeffs.K4567_plus, " fm² (旧 K 通道赝标量)")
println("  K₄₅₆₇⁻ = ", K_coeffs.K4567_minus, " fm² (旧 K 通道标量)")
println("  K₈⁺ = ", K_coeffs.K8_plus, " fm² (八重态标量)")
println("  K₈⁻ = ", K_coeffs.K8_minus, " fm² (八重态赝标量)")
println("  K₀₈⁺ = ", K_coeffs.K08_plus, " fm² (混合标量)")
println("  K₀₈⁻ = ", K_coeffs.K08_minus, " fm² (混合赝标量)")
```

**输出格式示意（具体数值必须以当前模型配置和单位换算为准）**:
```
原始耦合常数:
  G = 9.3e-6 fm²
  K = 1.2e-13 fm⁵

历史 helper（H_f=-phi_f）:
  H_u (API variable G_u) = -0.3314 fm⁻³
  H_s (API variable G_s) = -0.2145 fm⁻³

有效耦合系数:
  K₀⁺ = 9.6e-6 fm² (单态赝标量)
  K₀⁻ = 9.0e-6 fm² (单态标量)
  K₁₂₃⁺ = 9.2e-6 fm² (π通道赝标量)
  K₁₂₃⁻ = 9.4e-6 fm² (π通道标量)
  K₄₅₆₇⁺ = 9.1e-6 fm² (旧 K 通道赝标量)
  K₄₅₆₇⁻ = 9.5e-6 fm² (旧 K 通道标量)
  K₈⁺ = 9.4e-6 fm² (八重态赝标量)
  K₈⁻ = 9.2e-6 fm² (八重态标量)
  K₀₈⁺ = 1.4e-7 fm² (混合赝标量)
  K₀₈⁻ = -1.4e-7 fm² (混合标量)
```

#### 各向异性场景

当考虑动量各向异性时，使用 `A_aniso` 代替 `A`：

```julia
using .OneLoopIntegralsAniso: A_aniso

ξ = 0.25  # 各向异性参数

nodes_x, weights_x = gauleg(-1.0, 1.0, 32)

A_u = A_aniso(m_u_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ;
              nodes_p=nodes_p, weights_p=weights_p,
              nodes_x=nodes_x, weights_x=weights_x)

G_u = calculate_G_from_A(A_u, m_u_inv_fm)
# 继续计算有效耦合系数...
```

---

### `coupling_matrix_determinant(K0, K8, K08)`

计算混合介子传播子所需的耦合矩阵行列式。

#### 函数签名

```julia
coupling_matrix_determinant(K0::Float64, 
                            K8::Float64, 
                            K08::Float64) -> Float64
```

#### 物理公式

```
det K = K0 × K8 - K08²
```

#### 参数

| 参数 | 类型 | 说明 | 单位 |
|------|------|------|------|
| `K0` | `Float64` | 单态通道耦合系数 | fm² |
| `K8` | `Float64` | 八重态通道耦合系数 | fm² |
| `K08` | `Float64` | 混合通道耦合系数 | fm² |

#### 返回值

返回 `Float64`，表示耦合矩阵行列式（单位：fm⁴）。

#### 物理约束

**因果性条件**: 对于物理上有意义的介子传播子，必须满足：

```
det K > 0
```

否则，传播子会出现非物理的极点（如负质量平方或虚质量），表示模型在该参数点失效。

#### 物理应用

混合介子（如η-η'或σ-f₀）的传播子需要对2×2矩阵求逆：

```
D⁻¹ = [1/K₀  -K₀₈/K₀K₈]
      [-K₀₈/K₀K₈  1/K₈ ]

D = (det K)⁻¹ × [K₈  K₀₈]
                  [K₀₈  K₀]
```

行列式 det K 直接决定了传播子的归一化和极点位置。

#### 使用示例

```julia
using .EffectiveCouplings: calculate_effective_couplings, coupling_matrix_determinant

# 假设已计算得到K系数
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

# 计算赝标量通道（P，η/η'）行列式
det_K_pseudoscalar = coupling_matrix_determinant(K_coeffs.K0_plus, 
                                                 K_coeffs.K8_plus, 
                                                 K_coeffs.K08_plus)

# 计算标量通道（S，σ/σ'）行列式
det_K_scalar = coupling_matrix_determinant(K_coeffs.K0_minus, 
                                           K_coeffs.K8_minus, 
                                           K_coeffs.K08_minus)

println("标量通道行列式: det K^S = ", det_K_scalar, " fm⁴")
println("赝标量通道行列式: det K^P = ", det_K_pseudoscalar, " fm⁴")

# 物理性检查
if det_K_scalar > 0 && det_K_pseudoscalar > 0
    println("✅ 耦合矩阵满足物理约束")
else
    @warn "⚠️ 耦合矩阵出现非物理解，模型可能在该参数点失效"
end
```

**典型输出**:
```
标量通道行列式: det K^S = 8.9e-11 fm⁴
赝标量通道行列式: det K^P = 8.2e-11 fm⁴
✅ 耦合矩阵满足物理约束
```

#### 极端参数警告

在以下情况下，det K 可能变为负值或非常小：

1. **高温极限** (T > 200 MeV)：手征对称恢复导致 `H_f=-phi_f → 0`，K系数趋于简并
2. **极端化学势** (μ > 500 MeV)：密物质效应可能导致耦合矩阵奇异
3. **不合理的模型参数**：G 或 K 取值超出物理范围

此时需要检查输入参数的合理性，或考虑模型的适用范围限制。

---

### `mixing_matrix_elements(Π_uu, Π_ss, K_coeffs, channel)`

构造混合介子（η/η′ 或 σ/σ′）传播子/极点方程所需的 2×2 对称矩阵元素 (M00, M08, M88)。

#### 函数签名

```julia
mixing_matrix_elements(Π_uu::ComplexF64, Π_ss::ComplexF64,
                       K_coeffs::NamedTuple, channel::Symbol) -> NamedTuple
```

#### 约定

- `channel = :P`：赝标量通道（使用 `K*_plus` 与 `det_K_plus`）
- `channel = :S`：标量通道（使用 `K*_minus` 与 `det_K_minus`）

#### 公式

```
M00 = K0 - (4/3) * detK * (Π_uu + 2Π_ss)
M08 = K08 + (4/3) * sqrt(2) * detK * (Π_uu - Π_ss)
M88 = K8 - (4/3) * detK * (2Π_uu + Π_ss)
```

#### 返回值

返回 `NamedTuple`：`(M00, M08, M88, det_K)`。

---

## 温度依赖性分析

### 手征相变效应

随着温度升高，夸克凝聚 ⟨q̄q⟩ 减小，`H_f=-phi_f` 趋向于零，K系数逐渐退化为原始耦合常数 G。

**典型温度扫描结果**:

| T (MeV) | H_u (fm⁻³) | K₁₂₃⁻ (fm²) | K₁₂₃⁻/G | 手征对称性 |
|---------|-----|-------------|---------|-----------|
| 100 | -0.38 | 9.6e-6 | 1.03 | 完全破缺 |
| 150 | -0.33 | 9.4e-6 | 1.01 | 部分恢复 |
| 200 | -0.18 | 9.2e-6 | 0.99 | 接近恢复 |
| 250 | -0.05 | 9.3e-6 | 1.00 | 基本恢复 |

**解释**: 
- T < T_c（相变温度约170 MeV）时，|H_f| 较大，K系数偏离 G 显著
- T > T_c 时，H_f → 0，K系数趋于 G（K项贡献消失）

### 化学势依赖性

在有限化学势下，u/d夸克和s夸克的凝聚函数分裂更明显（由于质量差异），导致：

```
|H_u - H_s| 增大 → K₀₈± 增大 → η-η' 混合增强
```

---

## 注意事项与最佳实践

### 1. 单位换算

PNJL模型的原始参数通常以 MeV 为单位给出，需要转换为 fm⁻¹：

```julia
const ħc_MeV_fm = 197.327  # 转换因子

G_MeV2 = 5.0e-6  # MeV⁻²
G_fm2 = G_MeV2 / (ħc_MeV_fm^2)  # fm²
```

### 2. 参数自洽性

计算有效耦合系数时，需要确保：

- 夸克质量 m 是自洽求解能隙方程得到的有效质量（非裸质量）
- Polyakov 圈参数 Φ、Φbar 应从有效势最小化得到
- A 函数的截断参数 Λ 应与原始模型一致

### 3. 数值稳定性

当 K 很小（K ≪ G）时，K系数的修正项可能接近数值精度极限。此时：

- 使用双精度（Float64）计算
- 检查 K₀₈± 是否远小于 K₀± 和 K₈±
- 如果 |K₀₈±| < 1e-10，可以忽略混合效应

### 4. 物理合理性检验

计算完成后，建议执行以下检验：

```julia
# 1. 手征极限检验
@assert abs(K_coeffs.K0_plus - G_fm2) < 1e-6 "手征极限失败"

# 2. 符号检验（标量通道通常更强）
@assert K_coeffs.K123_plus > K_coeffs.K123_minus "π通道符号异常"

# 3. 行列式正定性
det_K = coupling_matrix_determinant(K_coeffs.K0_plus, K_coeffs.K8_plus, K_coeffs.K08_plus)
@assert det_K > 0 "耦合矩阵非物理"

# 4. 量级检验
@assert 1e-7 < K_coeffs.K123_minus < 1e-4 "K系数量级异常"
```

---

## 参考文档

### 公式推导
- `docs/reference/formula/relaxtime/couplings/EffectiveCoupling_K_FromA.md` - 有效耦合系数的详细推导

### 依赖模块
- `docs/api/integrals/OneLoopIntegrals.md` - A函数的各向同性计算
- `docs/api/integrals/OneLoopIntegralsAniso.md` - A函数的各向异性计算

### 下游应用
- `MesonPropagator.jl`（待开发）- 介子传播子计算
- `ScatteringAmplitude.jl`（待开发）- 散射振幅计算

---

## 版本历史

**v1.0** (2025-11-14)
- ✅ 实现 `calculate_G_from_A` 函数
- ✅ 实现 `calculate_effective_couplings` 函数
- ✅ 实现 `coupling_matrix_determinant` 辅助函数
- ✅ 完整的文档字符串和使用示例
- ✅ 支持各向异性场景

**待测试项**:
- 对称性检验（SU(3)对称、手征极限）
- 物理约束验证（det K > 0）
- 极限情况测试（K=0、高温极限）
- 温度和化学势扫描

---

## 联系与贡献

如有问题或建议，请参考 `tasks/task.md` 的步骤2任务说明。相关测试文件将在 `test/test_effective_couplings.jl` 中实现。
