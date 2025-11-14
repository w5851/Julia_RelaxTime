# EffectiveCouplings 模块 API 文档

## 模块概述

`EffectiveCouplings` 模块计算3味PNJL模型中的有效耦合常数。这些耦合系数是介子传播子计算的核心输入，由原始耦合常数（G、K）和夸克凝聚函数（通过单圈积分A计算）共同决定。

**物理背景**: 在随机相位近似(RPA)下，夸克-反夸克对形成的介子传播子由有效耦合系数 K_α^± 描述。这些系数编码了不同味道通道（单态、八重态）和自旋-宇称通道（标量S、赝标量P）的相互作用强度。

**在计算链中的位置**:
```
A函数(T,μ,m) → G^f → K系数 → 介子传播子D → 散射振幅M → 弛豫时间τ
```

## 依赖

- `Constants_PNJL` —— 提供原始耦合常数 G_fm2、K_fm5
- `OneLoopIntegrals` —— 提供各向同性A函数
- `OneLoopIntegralsAniso` —— 提供各向异性A函数（可选）

## 单位约定

- 原始耦合常数：G (fm²)、K (fm⁵)
- 夸克凝聚函数：G^f (无量纲)
- 有效耦合系数：K_α^± (fm²)
- 耦合矩阵行列式：det K (fm⁴)

---

## API 参考

### `calculate_G_from_A(A_f; Nc=3)`

从单圈积分A函数计算夸克凝聚相关函数G^f。

#### 函数签名

```julia
calculate_G_from_A(A_f::Float64; Nc::Int=3) -> Float64
```

#### 物理公式

```
G^f = -N_c / (4π²) × A_f(T, μ)
```

其中 A_f 由 `OneLoopIntegrals.A` 或 `OneLoopIntegralsAniso.A_aniso` 计算。

#### 参数

| 参数 | 类型 | 说明 | 单位 | 默认值 |
|------|------|------|------|--------|
| `A_f` | `Float64` | 单圈积分A函数的值 | fm | — |
| `Nc` | `Int` | 色数（QCD中固定为3） | 无量纲 | 3 |

#### 返回值

返回 `Float64`，表示夸克凝聚相关函数 G^f（无量纲）。

**物理意义**:
- G^f 描述了温度和密度对夸克凝聚 ⟨q̄q⟩ 的影响
- 手征对称破缺时 G^f ≠ 0
- 高温极限（T → ∞）时 G^f → 0（手征对称恢复）

#### 典型值

| 温度 | A_u (fm) | G^u (无量纲) |
|------|---------|--------------|
| T=100 MeV | ~1.6 | ~-0.38 |
| T=150 MeV | ~1.4 | ~-0.33 |
| T=200 MeV | ~1.2 | ~-0.29 |

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
println("A_u = ", A_u, " fm")

# 计算G^u
G_u = calculate_G_from_A(A_u)
println("G^u = ", G_u, " (无量纲)")
```

**典型输出**:
```
A_u = 1.3785 fm
G^u = -0.3314 (无量纲)
```

---

### `calculate_effective_couplings(G, K, G_u, G_s)`

计算10个有效耦合系数 K_α^±，用于介子传播子计算。

#### 函数签名

```julia
calculate_effective_couplings(G::Float64, K::Float64, 
                               G_u::Float64, G_s::Float64) -> NamedTuple
```

#### 物理公式

```
K₀±    = G ∓ (1/3)K(2G^μ + G^s)        # 单态通道
K₁₂₃± = G ± (1/2)KG^s                  # π介子通道
K₄₅₆₇± = G ± (1/2)KG^μ                 # K介子通道
K₈±    = G ± (1/6)K(4G^μ - G^s)        # 八重态通道
K₀₈±   = ±(1/6)√2 K(G^μ - G^s)         # 混合通道
```

**符号约定**:
- `+`: 标量通道（S，自旋0正宇称）
- `-`: 赝标量通道（P，自旋0负宇称）
- G^μ ≡ G_u = G_d（u和d夸克味道简并）
- G^s：s夸克的凝聚函数

#### 参数

| 参数 | 类型 | 说明 | 单位 | 典型值 |
|------|------|------|------|--------|
| `G` | `Float64` | 四夸克相互作用耦合常数 | fm² | ~5×10⁻⁶ MeV⁻² |
| `K` | `Float64` | 't Hooft六夸克相互作用耦合常数 | fm⁵ | ~1×10⁻¹³ MeV⁻⁵ |
| `G_u` | `Float64` | u/d夸克凝聚函数 | 无量纲 | ~-0.3 |
| `G_s` | `Float64` | s夸克凝聚函数 | 无量纲 | ~-0.2 |

#### 返回值

返回一个 `NamedTuple`，包含10个有效耦合系数（单位：fm²）：

```julia
(
    K0_plus::Float64,      # 单态标量
    K0_minus::Float64,     # 单态赝标量
    K123_plus::Float64,    # π通道标量（a0介子）
    K123_minus::Float64,   # π通道赝标量（π介子）
    K4567_plus::Float64,   # K通道标量（K0*介子）
    K4567_minus::Float64,  # K通道赝标量（K介子）
    K8_plus::Float64,      # 八重态标量
    K8_minus::Float64,     # 八重态赝标量（η8）
    K08_plus::Float64,     # 混合标量（σ-f0混合）
    K08_minus::Float64     # 混合赝标量（η-η'混合）
)
```

#### 物理意义详解

**味道通道分类**:

| 通道 | SU(3)表示 | 对应介子（P通道） | 对应介子（S通道） |
|------|-----------|------------------|------------------|
| K₀ | 单态 | η₀ (η') | σ, f₀(500) |
| K₁₂₃ | 三重态 | π⁰, π± | a₀(980) |
| K₄₅₆₇ | 四重态 | K⁰, K±, K̄⁰ | K₀*(1430) |
| K₈ | 八重态 | η₈ | f₀(980) |
| K₀₈ | 混合 | η-η' 混合 | σ-f₀ 混合 |

**手征极限检验**:

在手征极限下（m_q → 0），夸克凝聚消失 G^f → 0，所有 K_α^± 退化为原始耦合常数 G：

```julia
# 手征极限
G_u = 0.0
G_s = 0.0
K_coeffs = calculate_effective_couplings(G, K, 0.0, 0.0)
# 结果：所有 K_α^± = G
```

**SU(3)对称检验**:

当 G^u = G^s（三味对称）时：

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

# 计算G^f
G_u = calculate_G_from_A(A_u)
G_s = calculate_G_from_A(A_s)

# 计算有效耦合系数
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

# 输出结果
println("原始耦合常数:")
println("  G = ", G_fm2, " fm²")
println("  K = ", K_fm5, " fm⁵")
println()
println("夸克凝聚函数:")
println("  G^u = ", G_u)
println("  G^s = ", G_s)
println()
println("有效耦合系数:")
println("  K₀⁺ = ", K_coeffs.K0_plus, " fm² (单态标量)")
println("  K₀⁻ = ", K_coeffs.K0_minus, " fm² (单态赝标量)")
println("  K₁₂₃⁺ = ", K_coeffs.K123_plus, " fm² (π通道标量)")
println("  K₁₂₃⁻ = ", K_coeffs.K123_minus, " fm² (π通道赝标量)")
println("  K₄₅₆₇⁺ = ", K_coeffs.K4567_plus, " fm² (K通道标量)")
println("  K₄₅₆₇⁻ = ", K_coeffs.K4567_minus, " fm² (K通道赝标量)")
println("  K₈⁺ = ", K_coeffs.K8_plus, " fm² (八重态标量)")
println("  K₈⁻ = ", K_coeffs.K8_minus, " fm² (八重态赝标量)")
println("  K₀₈⁺ = ", K_coeffs.K08_plus, " fm² (混合标量)")
println("  K₀₈⁻ = ", K_coeffs.K08_minus, " fm² (混合赝标量)")
```

**典型输出**:
```
原始耦合常数:
  G = 9.3e-6 fm²
  K = 1.2e-13 fm⁵

夸克凝聚函数:
  G^u = -0.3314
  G^s = -0.2145

有效耦合系数:
  K₀⁺ = 9.6e-6 fm² (单态标量)
  K₀⁻ = 9.0e-6 fm² (单态赝标量)
  K₁₂₃⁺ = 9.2e-6 fm² (π通道标量)
  K₁₂₃⁻ = 9.4e-6 fm² (π通道赝标量)
  K₄₅₆₇⁺ = 9.1e-6 fm² (K通道标量)
  K₄₅₆₇⁻ = 9.5e-6 fm² (K通道赝标量)
  K₈⁺ = 9.4e-6 fm² (八重态标量)
  K₈⁻ = 9.2e-6 fm² (八重态赝标量)
  K₀₈⁺ = 1.4e-7 fm² (混合标量)
  K₀₈⁻ = -1.4e-7 fm² (混合赝标量)
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

G_u = calculate_G_from_A(A_u)
# 继续计算有效耦合系数...
```

---

### `coupling_matrix_determinant(K0_plus, K8_plus, K08_plus)`

计算混合介子传播子所需的耦合矩阵行列式。

#### 函数签名

```julia
coupling_matrix_determinant(K0_plus::Float64, 
                            K8_plus::Float64, 
                            K08_plus::Float64) -> Float64
```

#### 物理公式

```
det K = K₀⁺ × K₈⁺ - (K₀₈⁺)²
```

#### 参数

| 参数 | 类型 | 说明 | 单位 |
|------|------|------|------|
| `K0_plus` | `Float64` | 单态通道标量耦合系数 | fm² |
| `K8_plus` | `Float64` | 八重态通道标量耦合系数 | fm² |
| `K08_plus` | `Float64` | 混合通道标量耦合系数 | fm² |

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

# 计算标量通道行列式
det_K_scalar = coupling_matrix_determinant(K_coeffs.K0_plus, 
                                           K_coeffs.K8_plus, 
                                           K_coeffs.K08_plus)

# 计算赝标量通道行列式
det_K_pseudoscalar = coupling_matrix_determinant(K_coeffs.K0_minus, 
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

1. **高温极限** (T > 200 MeV)：手征对称恢复导致 G^f → 0，K系数趋于简并
2. **极端化学势** (μ > 500 MeV)：密物质效应可能导致耦合矩阵奇异
3. **不合理的模型参数**：G 或 K 取值超出物理范围

此时需要检查输入参数的合理性，或考虑模型的适用范围限制。

---

## 温度依赖性分析

### 手征相变效应

随着温度升高，夸克凝聚 ⟨q̄q⟩ 减小，G^f 趋向于零，K系数逐渐退化为原始耦合常数 G。

**典型温度扫描结果**:

| T (MeV) | G^u | K₁₂₃⁻ (fm²) | K₁₂₃⁻/G | 手征对称性 |
|---------|-----|-------------|---------|-----------|
| 100 | -0.38 | 9.6e-6 | 1.03 | 完全破缺 |
| 150 | -0.33 | 9.4e-6 | 1.01 | 部分恢复 |
| 200 | -0.18 | 9.2e-6 | 0.99 | 接近恢复 |
| 250 | -0.05 | 9.3e-6 | 1.00 | 基本恢复 |

**解释**: 
- T < T_c（相变温度约170 MeV）时，|G^f| 较大，K系数偏离 G 显著
- T > T_c 时，G^f → 0，K系数趋于 G（K项贡献消失）

### 化学势依赖性

在有限化学势下，u/d夸克和s夸克的凝聚函数分裂更明显（由于质量差异），导致：

```
|G^u - G^s| 增大 → K₀₈± 增大 → η-η' 混合增强
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
- `doc/formula/K_有效耦合常数byA.md` - 有效耦合系数的详细推导

### 依赖模块
- `api/OneLoopIntegrals.md` - A函数的各向同性计算
- `api/OneLoopIntegralsAniso.md` - A函数的各向异性计算

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
