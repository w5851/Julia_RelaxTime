# MesonPropagator 模块 API 文档

## 模块概述

`MesonPropagator` 模块计算3味PNJL模型中的介子传播子，支持一般介子（π、K、σ_π、σ_K）和混合介子（η/η'、σ/σ'）的完整计算。传播子是计算散射振幅和弛豫时间的核心物理量。

**物理背景**: 在随机相位近似(RPA)下，夸克-反夸克对通过单圈图形成介子束缚态。介子传播子描述了这些束缚态的传播行为，其极点位置对应介子质量，留数对应耦合强度。

**在计算链中的位置**:
```
K系数 → 极化函数Π → 介子传播子D → 散射振幅M → 弛豫时间τ
```

## 依赖

- `Constants_PNJL` —— 提供Gell-Mann矩阵（λ₀、λ₈）和夸克味波函数（ψ_u、ψ_d、ψ_s及其共轭）
- `EffectiveCouplings` —— 提供K系数计算和det_K函数
- `LinearAlgebra` —— 提供矩阵运算（det、inv、transpose）

## 单位约定

- 有效耦合系数K：fm²
- 极化函数Π：fm²（ComplexF64）
- 耦合矩阵M：无量纲（ComplexF64矩阵）
- det_K：fm⁴
- 传播子D：fm²（ComplexF64）

## 各向异性修正的设计哲学

**重要**: 本模块不直接依赖各向异性参数ξ，这是有意的模块化设计：

- 各向异性通过预计算的极化函数Π间接影响传播子
- 调用链：`PolarizationAniso.polarization_aniso(ξ, ...)` → Π(ComplexF64) → `meson_propagator(...)`
- 使用时在调用链最前端选择`Polarization.jl`(各向同性)或`PolarizationAniso.jl`(各向异性)
- 传播子函数无需修改，保持模块简洁性

---

## API 参考

### 主函数

#### `meson_propagator_simple(meson_type, K_coeffs, Π)`

计算一般介子（π、K、σ_π、σ_K）的传播子。

##### 函数签名

```julia
meson_propagator_simple(meson_type::Symbol, K_coeffs::NamedTuple, Π::ComplexF64) -> ComplexF64
```

##### 物理公式

```
D(q²) = 1 / (1 - K_α^± Π(q²))
```

其中α为介子通道，±表示自旋-宇称通道。

##### 介子类型映射表

| 介子 | 通道 | 夸克组合 | 极化函数 | K系数变量名 |
|------|------|----------|---------|-------------|
| π | P | ūu/d̄d | Π_{uu}^P | `K123_plus` |
| K | P | ūs/d̄s | Π_{us}^P | `K4567_plus` |
| σ_π | S | ūu/d̄d | Π_{uu}^S | `K123_minus` |
| σ_K | S | ūs/d̄s | Π_{us}^S | `K4567_minus` |

##### K系数命名对照表

| 物理符号 | 代码变量名 | 说明 | 对应介子 |
|---------|-----------|------|----------|
| K₀± | `K0_plus/minus` | 单态通道 | η'/σ |
| K₁±=K₂±=K₃± | `K123_plus/minus` | π通道(三个分量相等) | π, a₀ |
| K₄±=K₅±=K₆±=K₇± | `K4567_plus/minus` | K通道(四个分量相等) | K, K₀* |
| K₈± | `K8_plus/minus` | 八重态通道 | η₈, f₀ |
| K₀₈± | `K08_plus/minus` | 混合通道 | η-η'混合 |

**关键映射规则**:
- 赝标量P（π、K）使用K⁺系数（`*_plus`）
- 标量S（σ_π、σ_K）使用K⁻系数（`*_minus`）

##### K系数自动选择逻辑

函数内部根据`meson_type`自动选择正确的K系数，降低调用者使用难度：

| meson_type | 自动选择的K系数 | 物理原因 |
|------------|----------------|----------|
| `:pi` | `K123_plus` | π是赝标量P通道，用K⁺ |
| `:K` | `K4567_plus` | K是赝标量P通道，用K⁺ |
| `:sigma_pi` | `K123_minus` | σ_π是标量S通道，用K⁻ |
| `:sigma_K` | `K4567_minus` | σ_K是标量S通道，用K⁻ |

##### 参数

| 参数 | 类型 | 说明 | 单位 |
|------|------|------|------|
| `meson_type` | `Symbol` | 介子类型（`:pi`, `:K`, `:sigma_pi`, `:sigma_K`） | — |
| `K_coeffs` | `NamedTuple` | 预计算的K系数（通过`EffectiveCouplings.calculate_effective_couplings`获取） | fm² |
| `Π` | `ComplexF64` | 预计算的极化函数（注意：K介子使用Π_{us}，不是Π_{uu}） | fm² |

##### 返回值

返回 `ComplexF64`，表示介子传播子D（单位：fm²）。

**物理意义**:
- 实部：传播子的色散关系
- 虚部：介子的衰变宽度（若存在）
- 极点位置（Re D → ∞）：介子质量
- 留数：耦合强度

##### 典型值

| 介子 | q² (fm⁻²) | |D| (fm²) | 物理状态 |
|------|-----------|-----------|----------|
| π | -0.02 | ~100 | 接近质量壳 |
| π | -0.1 | ~10 | 远离质量壳 |
| K | -0.025 | ~80 | 接近质量壳 |

##### 使用示例

```julia
using .Constants_PNJL: G_fm2, K_fm5, ħc_MeV_fm
using .EffectiveCouplings: calculate_effective_couplings
using .MesonPropagator: meson_propagator_simple

# 预计算K系数（所有介子共用，批量计算时只需计算一次）
G_u = -0.3  # 通过calculate_G_from_A得到
G_s = -0.2
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

# 预计算极化函数（通过Polarization或PolarizationAniso模块）
Π_uu_P = 1.0e-5 + 0.0im  # ūu极化函数（赝标量通道）
Π_us_P = 9.0e-6 + 0.0im  # ūs极化函数（赝标量通道）
Π_uu_S = 1.2e-5 + 0.0im  # ūu极化函数（标量通道）
Π_us_S = 1.0e-5 + 0.0im  # ūs极化函数（标量通道）

# 计算不同介子的传播子
D_pi = meson_propagator_simple(:pi, K_coeffs, Π_uu_P)
D_K = meson_propagator_simple(:K, K_coeffs, Π_us_P)  # 注意使用Π_us
D_sigma_pi = meson_propagator_simple(:sigma_pi, K_coeffs, Π_uu_S)
D_sigma_K = meson_propagator_simple(:sigma_K, K_coeffs, Π_us_S)

println("D_π = ", D_pi, " fm²")
println("D_K = ", D_K, " fm²")
println("D_σπ = ", D_sigma_pi, " fm²")
println("D_σK = ", D_sigma_K, " fm²")
```

---

#### `meson_propagator_mixed(det_K, M_matrix, q1, q2, q3, q4, channel)`

计算混合介子（η/η'、σ/σ'）的传播子。

##### 函数签名

```julia
meson_propagator_mixed(det_K::Float64, M_matrix::Matrix{ComplexF64},
                      q1::Symbol, q2::Symbol, q3::Symbol, q4::Symbol,
                      channel::Symbol) -> ComplexF64
```

##### 物理公式

```
D(q²) = [2det(K) / det(M)] × J^T M J'
```

其中：
- det(K)：耦合矩阵行列式（预计算）
- M：2×2耦合矩阵（包含极化函数修正）
- J, J'：流算符向量（由散射过程的味量子数决定）

##### 散射道映射表

| 散射道 | 物理过程 | ψ | ψ̄ | ψ' | ψ̄' | 介子四动量 |
|-------|---------|---|---|----|----|----------|
| t道 | q₁+q₂→q₃+q₄(纯夸克) | q1 | q3 | q2 | q4 | q=p₁-p₃ |
| s道 | q₁+q̄₂→q₃+q̄₄(湮灭) | q1 | q2 | q3 | q4 | P=p₁+p₂ |
| u道 | q₁+q₂→q₃+q₄(交叉) | q1 | q4 | q2 | q3 | q=p₁-p₄ |

**场算符的产生/湮灭性质**（关键物理概念）：
- 对于**夸克**：ψ是**湮灭算符**，ψ̄是**产生算符**
- 对于**反夸克**：ψ是**产生算符**，ψ̄是**湮灭算符**
- 判断依据：费曼图中箭头指向顶点为湮灭，离开顶点为产生

**散射过程约束**：
- q1和q3必定是夸克（`:u`, `:d`, `:s`）
- q2和q4要么全为夸克，要么全为反夸克
- t道和u道：纯夸克散射
- s道：夸克-反夸克湮灭

##### 参数

| 参数 | 类型 | 说明 | 单位 |
|------|------|------|------|
| `det_K` | `Float64` | 预计算的耦合矩阵行列式（通过`EffectiveCouplings.coupling_matrix_determinant`获取） | fm⁴ |
| `M_matrix` | `Matrix{ComplexF64}` | 预计算的2×2复数耦合矩阵（通过`calculate_coupling_matrix`获取） | 无量纲 |
| `q1, q2, q3, q4` | `Symbol` | 散射过程q₁+q₂→q₃+q₄的夸克/反夸克类型（`:u`, `:d`, `:s`, `:ubar`, `:dbar`, `:sbar`） | — |
| `channel` | `Symbol` | 散射道类型（`:t`, `:s`, `:u`） | — |

##### 返回值

返回 `ComplexF64`，表示混合介子传播子D（单位：fm²）。

##### 使用示例

```julia
using .Constants_PNJL: G_fm2, K_fm5
using .EffectiveCouplings: calculate_effective_couplings, coupling_matrix_determinant
using .MesonPropagator: meson_propagator_mixed, calculate_coupling_matrix

# 预计算K系数
G_u = -0.3
G_s = -0.2
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

# 预计算det_K（赝标量通道，用于η/η'）
det_K_P = coupling_matrix_determinant(K_coeffs.K0_plus, K_coeffs.K8_plus, K_coeffs.K08_plus)

# 预计算极化函数和M矩阵
Π_uu_P = 1.0e-5 + 0.0im
Π_ss_P = 8.0e-6 + 0.0im
M_P = calculate_coupling_matrix(Π_uu_P, Π_ss_P, K_coeffs, :P)

# 计算η/η'传播子（u+d→u+d散射，t道）
D_eta = meson_propagator_mixed(det_K_P, M_P, :u, :d, :u, :d, :t)
println("D_η/η' (t道) = ", D_eta, " fm²")

# 计算η/η'传播子（u+d̅→u+d̅散射，s道）
D_eta_s = meson_propagator_mixed(det_K_P, M_P, :u, :dbar, :u, :dbar, :s)
println("D_η/η' (s道) = ", D_eta_s, " fm²")

# 计算σ/σ'传播子（标量通道）
det_K_S = coupling_matrix_determinant(K_coeffs.K0_minus, K_coeffs.K8_minus, K_coeffs.K08_minus)
M_S = calculate_coupling_matrix(Π_uu_P, Π_ss_P, K_coeffs, :S)
D_sigma = meson_propagator_mixed(det_K_S, M_S, :u, :d, :u, :d, :t)
println("D_σ/σ' (t道) = ", D_sigma, " fm²")
```

---

### 辅助函数

#### `calculate_coupling_matrix(Π_uu, Π_ss, K_coeffs, channel)`

计算混合介子传播子所需的2×2复数耦合矩阵M。

##### 函数签名

```julia
calculate_coupling_matrix(Π_uu::ComplexF64, Π_ss::ComplexF64, 
                         K_coeffs::NamedTuple, channel::Symbol) -> Matrix{ComplexF64}
```

##### 物理公式

M矩阵元素：
```
M₀₀ = 1 - K₀± Π_{uu}
M₀₈ = M₈₀ = -K₀₈± Π_{uu} × 4√2/3
M₈₈ = 1 - K₈± [Π_{uu} × 4/3 + Π_{ss} × 2/3]
```

其中±表示通道类型（+为赝标量P，-为标量S）。

##### 参数

| 参数 | 类型 | 说明 | 单位 |
|------|------|------|------|
| `Π_uu` | `ComplexF64` | u/d夸克的极化函数 | fm² |
| `Π_ss` | `ComplexF64` | s夸克的极化函数 | fm² |
| `K_coeffs` | `NamedTuple` | 预计算的K系数 | fm² |
| `channel` | `Symbol` | 通道类型（`:P`为赝标量，`:S`为标量） | — |

##### 返回值

返回2×2复数矩阵M（无量纲）。

##### 重要修正

M₀₈系数使用`4*sqrt(2)/3` ≈ 1.8856，**不是**4/(3√2) ≈ 0.9428。这是根据公式文档的修正。

##### 使用示例

```julia
# 计算赝标量通道（η/η'）的M矩阵
M_P = calculate_coupling_matrix(Π_uu_P, Π_ss_P, K_coeffs, :P)

# 计算标量通道（σ/σ'）的M矩阵
M_S = calculate_coupling_matrix(Π_uu_S, Π_ss_S, K_coeffs, :S)
```

---

#### `extract_flavor(q)` 和 `get_quark_wavefunction(flavor, is_bar)`

从夸克/反夸克符号中提取味类型并返回对应的波函数（内部辅助函数）。

##### 函数签名

```julia
extract_flavor(q::Symbol) -> (Symbol, Bool)
get_quark_wavefunction(flavor::Symbol, is_bar::Bool) -> Union{Vector{Float64}, Matrix{Float64}}
```

##### 参数与返回值

- `extract_flavor(:ubar)` 返回 `(:u, true)`（u味的反夸克）
- `get_quark_wavefunction(:u, false)` 返回列向量 `[1.0, 0.0, 0.0]`
- `get_quark_wavefunction(:u, true)` 返回行向量 `[1.0 0.0 0.0]`（1×3矩阵）

---

#### `calculate_current_vector(q1, q2, channel)`

根据散射道映射表自动选择场算符，计算流算符向量（内部辅助函数）。

##### 函数签名

```julia
calculate_current_vector(q1::Symbol, q2::Symbol, channel::Symbol) -> Vector{Float64}
```

##### 物理公式

```
J = [ψ̄·λ₀·ψ, ψ̄·λ₈·ψ]
```

其中矩阵乘法顺序：ψbar (1×3) × λ (3×3) × ψ (3×1) = 标量。

##### 返回值

返回2×1列向量 J = [J₀, J₈]。

**物理意义**：
- J₀：单态分量
- J₈：八重态分量
- 流算符向量编码了介子的味量子数

---

## Gell-Mann矩阵和味波函数常量

本模块使用`Constants_PNJL`中预定义的常量：

### Gell-Mann矩阵（3×3矩阵）

```julia
λ₀  # 味单位矩阵（归一化）
λ₈  # Gell-Mann第8矩阵（u,d对称，s不同）
```

### 夸克味波函数

**列向量**（夸克）：
```julia
ψ_u = [1.0, 0.0, 0.0]
ψ_d = [0.0, 1.0, 0.0]
ψ_s = [0.0, 0.0, 1.0]
```

**行向量**（反夸克，1×3矩阵）：
```julia
ψbar_u = [1.0 0.0 0.0]
ψbar_d = [0.0 1.0 0.0]
ψbar_s = [0.0 0.0 1.0]
```

**命名约定**：代码中使用ASCII命名（如`ψbar_u`），避免Unicode兼容性问题。

---

## 旋量处理说明

**重要概念**：本模块只处理味空间的3×3矩阵运算。

- Dirac旋量指标已通过矩阵乘法$\bar{\psi}\lambda\psi$自动收缩
- 公式中的ψ和ψ̄表示夸克场，但在计算中只保留味空间部分
- 矩阵乘法顺序：从左到右依次进行（1×3矩阵 × 3×3矩阵 × 3×1向量）

---

## 性能优化建议

### 批量计算策略

当计算多个介子或多个动量点的传播子时，应采用以下优化策略：

1. **预计算K系数**（所有介子共用）
```julia
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
```

2. **预计算det_K和M矩阵**（混合介子共用）
```julia
det_K_P = coupling_matrix_determinant(K_coeffs.K0_plus, K_coeffs.K8_plus, K_coeffs.K08_plus)
M_P = calculate_coupling_matrix(Π_uu_P, Π_ss_P, K_coeffs, :P)
```

3. **循环中只计算极化函数和传播子**
```julia
for q2 in q2_values
    Π_uu = polarization_function(q2, ...)  # 各向同性或各向异性
    D_pi = meson_propagator_simple(:pi, K_coeffs, Π_uu)
    # ... 其他计算
end
```

**预期性能提升**：3-5倍（相比每次重新计算K系数）

---

## 物理约束与注意事项

### 物理合理性检验

1. **质量壳条件**：q² → -m²时，|D| → ∞（传播子极点）
2. **无衰变条件**：若介子稳定，Im D ≈ 0
3. **因果性**：det K > 0（保证传播子因果性）
4. **手征极限**：G^μ = G^s = 0时，所有K_α^±退化为G

### 常见错误

1. **K介子极化函数错误**：必须使用Π_{us}，不是Π_{uu}
2. **通道-K系数映射错误**：赝标量P用K⁺，标量S用K⁻（与直觉相反）
3. **M₀₈系数错误**：使用`4*sqrt(2)/3`，不是`4/(3√2)`
4. **散射道映射错误**：未正确根据channel选择场算符对应关系

---

## 与其他模块的集成

### 完整计算流程

```julia
# 第1步：计算A函数（OneLoopIntegrals或OneLoopIntegralsAniso）
A_u = A(m_u, μ, T, Φ, Φbar, nodes, weights)
A_s = A(m_s, μ, T, Φ, Φbar, nodes, weights)

# 第2步：计算G^f（EffectiveCouplings）
G_u = calculate_G_from_A(A_u)
G_s = calculate_G_from_A(A_s)

# 第3步：计算K系数（EffectiveCouplings）
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

# 第4步：计算极化函数（Polarization或PolarizationAniso）
Π_uu_P = polarization_uu(q2, :P, ...)  # 各向同性
# 或
Π_uu_P = polarization_aniso(ξ, q2, :P, ...)  # 各向异性

# 第5步：计算传播子（MesonPropagator）
D_pi = meson_propagator_simple(:pi, K_coeffs, Π_uu_P)
```

### 各向异性计算示例

```julia
# 各向同性（ξ=0）
Π_iso = polarization_uu(q2, :P, ...)
D_iso = meson_propagator_simple(:pi, K_coeffs, Π_iso)

# 各向异性（ξ≠0）
Π_aniso = polarization_aniso(ξ, q2, :P, ...)
D_aniso = meson_propagator_simple(:pi, K_coeffs, Π_aniso)

# 比较各向异性修正的影响
Δ_aniso = abs(D_aniso - D_iso) / abs(D_iso)
println("各向异性修正幅度: ", Δ_aniso * 100, "%")
```

---

## 已知限制

1. **衰变宽度计算**：当前版本未实现介子衰变宽度的自动提取
2. **质量提取**：未实现从传播子极点自动提取介子质量的功能
3. **Mott相变检测**：未实现det K符号变化的自动监测
4. **高阶修正**：仅包含单圈RPA图，未考虑高阶量子修正

这些功能留待后续步骤或作为独立工具函数实现。
