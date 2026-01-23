# OneLoopIntegralsAniso 模块 API 文档

## 模块概述

`OneLoopIntegralsAniso` 模块（正式名称为 `OneLoopIntegralsCorrection`）提供动量各向异性情况下的单圈积分修正计算。模块包含两类函数：
1. **一阶修正函数**：`A_correction` 和 `B0_correction`，计算各向异性参数 ξ 的一阶修正项
2. **完整各向异性函数**：`A_aniso`，直接计算包含 Romatschke-Strickland (RS) 分布的完整积分

本模块实现依据 `docs/reference/formula/relaxtime/integrals/OneLoopIntegral_A_Aniso_Collinear.md`、`docs/reference/formula/relaxtime/integrals/OneLoopIntegral_B0_Aniso_General.md` 等公式文档，适用于重离子碰撞等非平衡系统的早期动力学研究。

## 依赖

- `QuadGK.jl` —— 自适应一维积分器
- `OneLoopIntegrals` —— 各向同性基础积分（内部依赖）
- `PNJLQuarkDistributions_Aniso` —— 各向异性夸克分布函数
- `GaussLegendre` —— 高斯-勒让德节点生成

## 单位约定

遵循项目统一的自然单位制 (ℏ = c = 1)：

- 所有动力学参数（`m`, `μ`, `T`, `λ`, `k`）单位为 fm⁻¹
- 各向异性参数 `ξ` 为无量纲（典型取值范围 0~0.5）
- `A_correction` 和 `A_aniso` 返回值单位：fm
- `B0_correction` 返回值 `(real_part, imag_part)` 单位：fm

---

## 函数使用指南

### 一阶修正 vs 完整计算

**何时使用一阶修正**：
- ξ 较小（通常 ξ < 0.3）时，一阶修正可提供足够精度
- 需要快速计算，对性能要求较高
- 已有各向同性结果 `A`，只需叠加修正项：`A_total = A + A_correction`

**何时使用完整计算**：
- ξ 较大（ξ ≥ 0.3）时，一阶近似可能不足
- 需要最高精度的各向异性结果
- 直接计算完整 RS 分布下的积分：`A_total = A_aniso`

**性能对比**：
- `A + A_correction`：计算两次独立积分，但 `A` 可复用于多个 ξ 值
- `A_aniso`：单次计算，但需要双重积分（动量 + 角度），耗时约为前者 1.5-2 倍

---

## API 参考

### `A_correction(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p)`

计算单传播子积分 A 在动量各向异性下的**一阶修正项**。此函数返回的是 ξ 的线性响应，需与各向同性项 `A(m, μ, T, Φ, Φbar, nodes_p, weights_p)` 相加才能得到完整结果。

#### 函数签名

```julia
A_correction(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64,
             ξ::Float64, nodes_p::Vector{Float64}, weights_p::Vector{Float64}) -> Float64
```

#### 参数

| 参数 | 类型 | 说明 | 单位 | 典型值 |
|------|------|------|------|--------|
| `m` | `Float64` | 夸克有效质量 | fm⁻¹ | 0.2-0.4 |
| `μ` | `Float64` | 夸克化学势 | fm⁻¹ | 0-0.3 |
| `T` | `Float64` | 温度 | fm⁻¹ | 0.1-0.3 |
| `Φ` | `Float64` | Polyakov 圈参数 | 无量纲 | 0-1 |
| `Φbar` | `Float64` | 共轭 Polyakov 圈参数 | 无量纲 | 0-1 |
| `ξ` | `Float64` | 各向异性参数（0 表示各向同性） | 无量纲 | 0-0.5 |
| `nodes_p` | `Vector{Float64}` | 动量积分节点 | fm⁻¹ | 由 `gauleg(0, 20, 64)` 生成 |
| `weights_p` | `Vector{Float64}` | 动量积分权重 | fm⁻¹ | 与 `nodes_p` 配套 |

#### 返回值

返回 `Float64`，即 A 的一阶修正项 ΔA₁(ξ)。

**重要**：此值需与各向同性项相加：
```julia
A_isotropic = OneLoopIntegrals.A(m, μ, T, Φ, Φbar, nodes_p, weights_p)
A_total = A_isotropic + A_correction(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p)
```

#### 物理意义

在动量各向异性下，夸克分布函数修正为（一阶展开）：
```math
f^{\text{aniso}} \approx f^0 - \frac{\xi (p\cos\theta)^2}{2ET} \frac{df^0}{dE}
```

`A_correction` 计算上式第二项对 A 积分的贡献：
```math
\Delta A_1 = \frac{4}{3} \int dp \, p^2 \left[ -\frac{\xi p^2}{2ET} \left( \frac{df_q}{dE} + \frac{df_{\bar{q}}}{dE} \right) \right]
```

#### 数值实现要点

- 使用高斯-勒让德积分对动量 p 求积
- 角度项 `(cosθ)²` 的平均值为 1/3（已解析处理）
- ξ=0 时应严格返回 0（可用于验证实现正确性）
- 线性性质：`A_correction(2ξ) / A_correction(ξ) ≈ 2.0`（小 ξ 时）

#### 使用示例

```julia
using Pkg
Pkg.activate(".")

include("src/relaxtime/OneLoopIntegrals.jl")
include("src/relaxtime/OneLoopIntegralsAniso.jl")
include("src/integration/GaussLegendre.jl")

using .OneLoopIntegrals: A
using .OneLoopIntegralsCorrection: A_correction
using .GaussLegendre: gauleg

# 生成积分节点
nodes_p, weights_p = gauleg(0.0, 20.0, 64)

# 物理参数
m = 0.3      # fm⁻¹
μ = 0.1      # fm⁻¹
T = 0.15     # fm⁻¹
Φ = 0.2
Φbar = 0.2
ξ = 0.1      # 小各向异性

# 计算各向同性项和修正项
A_iso = A(m, μ, T, Φ, Φbar, nodes_p, weights_p)
ΔA = A_correction(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p)
A_total = A_iso + ΔA

println("A (各向同性): ", A_iso, " fm")
println("ΔA (一阶修正): ", ΔA, " fm")
println("A (总计):      ", A_total, " fm")
println("修正比例:      ", 100*ΔA/A_iso, " %")
```

**典型输出**（ξ=0.1 时）：
```
A (各向同性): 1.378 fm
ΔA (一阶修正): -0.032 fm
A (总计):      1.346 fm
修正比例:      -2.3 %
```

---

### `A_aniso(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p, nodes_cosθ, weights_cosθ)`

计算单传播子积分 A 在动量各向异性下的**完整形式**，使用 Romatschke-Strickland 分布函数直接积分，无需一阶近似。

#### 函数签名

```julia
A_aniso(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64,
        ξ::Float64, nodes_p::Vector{Float64}, weights_p::Vector{Float64},
        nodes_cosθ::Vector{Float64}, weights_cosθ::Vector{Float64}) -> Float64
```

#### 参数

| 参数 | 类型 | 说明 | 单位 | 典型值 |
|------|------|------|------|--------|
| `m` | `Float64` | 夸克有效质量 | fm⁻¹ | 0.2-0.4 |
| `μ` | `Float64` | 夸克化学势 | fm⁻¹ | 0-0.3 |
| `T` | `Float64` | 温度 | fm⁻¹ | 0.1-0.3 |
| `Φ` | `Float64` | Polyakov 圈参数 | 无量纲 | 0-1 |
| `Φbar` | `Float64` | 共轭 Polyakov 圈参数 | 无量纲 | 0-1 |
| `ξ` | `Float64` | 各向异性参数 | 无量纲 | 0-1 |
| `nodes_p` | `Vector{Float64}` | 动量积分节点 | fm⁻¹ | 由 `gauleg(0, 20, 64)` 生成 |
| `weights_p` | `Vector{Float64}` | 动量积分权重 | fm⁻¹ | - |
| `nodes_cosθ` | `Vector{Float64}` | 角度积分节点（cosθ ∈ [-1, 1]） | 无量纲 | 由 `gauleg(-1, 1, 32)` 生成 |
| `weights_cosθ` | `Vector{Float64}` | 角度积分权重 | 无量纲 | - |

#### 返回值

返回 `Float64`，即完整各向异性修正下的 A 积分值，可直接使用，无需额外叠加各向同性项。

#### 物理意义

使用完整的 RS 分布函数：
```math
f^{\text{RS}}(p, \theta) = f^0\left( \sqrt{p^2 + m^2 + \xi (p\cos\theta)^2} \right)
```

计算：
```math
A_{\text{aniso}} = 2 \int_0^\infty dp \int_{-1}^{1} d(\cos\theta) \, p^2 \frac{1}{E_{\text{aniso}}} \left[ f_q^{\text{RS}} + f_{\bar{q}}^{\text{RS}} \right] - A_{\text{const}}
```

其中：
- $E_{\text{aniso}} = \sqrt{p^2 + m^2 + \xi (p\cos\theta)^2}$
- $A_{\text{const}}$ 为常数项（与 `OneLoopIntegrals.A` 相同）
- 因子 2 来自对 ϕ 角积分（各向异性仅沿 z 轴）

#### 数值实现要点

- 双重积分：外层对动量 p，内层对 cosθ
- ξ=0 时应退化为各向同性结果：`A_aniso(ξ=0) ≈ A(isotropic)`
- 对于大 ξ（>0.5），此函数比一阶近似更可靠
- 角度积分节点数建议 ≥32（确保角度依赖充分捕捉）

#### 使用示例

```julia
using Pkg
Pkg.activate(".")

include("src/relaxtime/OneLoopIntegralsAniso.jl")
include("src/integration/GaussLegendre.jl")

using .OneLoopIntegralsCorrection: A_aniso
using .GaussLegendre: gauleg

# 生成积分节点（动量 + 角度）
nodes_p, weights_p = gauleg(0.0, 20.0, 64)
nodes_cosθ, weights_cosθ = gauleg(-1.0, 1.0, 32)

# 物理参数
m = 0.3
μ = 0.1
T = 0.15
Φ = 0.2
Φbar = 0.2
ξ = 0.4      # 较大各向异性

# 直接计算完整各向异性结果
A_full = A_aniso(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p, nodes_cosθ, weights_cosθ)

println("A (完整各向异性, ξ=", ξ, "): ", A_full, " fm")
```

**典型输出**（ξ=0.4 时）：
```
A (完整各向异性, ξ=0.4): 1.256 fm
```

#### 与 A_correction 的对比

```julia
# 方法 1：一阶近似
A_iso = A(m, μ, T, Φ, Φbar, nodes_p, weights_p)
ΔA = A_correction(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p)
A_approx = A_iso + ΔA

# 方法 2：完整计算
A_exact = A_aniso(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p, nodes_cosθ, weights_cosθ)

# 相对误差
rel_error = abs(A_exact - A_approx) / A_exact * 100
println("相对误差: ", rel_error, " %")
```

**预期结果**：
- ξ=0.1：相对误差 < 0.5%（一阶近似优秀）
- ξ=0.3：相对误差 ≈ 2-3%（一阶近似可用）
- ξ=0.5：相对误差 ≈ 8-10%（建议用完整计算）

---

### `B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ; rtol=1e-3, atol=0.0)`

计算双传播子积分 B₀ 在动量各向异性下的**一阶修正项**。返回修正后的实部和虚部，需与各向同性的 `B0` 结果相加。

#### 函数签名

```julia
B0_correction(λ::Float64, k::Float64, m1::Float64, m2::Float64,
              μ1::Float64, μ2::Float64, T::Float64,
              Φ::Float64, Φbar::Float64, ξ::Float64;
              rtol::Float64=1e-3, atol::Float64=0.0) -> NTuple{2, Float64}
```

#### 参数

| 参数 | 类型 | 说明 | 单位 | 默认值 |
|------|------|------|------|--------|
| `λ` | `Float64` | 能量-化学势组合参数 | fm⁻¹ | - |
| `k` | `Float64` | 外部三动量模长 | fm⁻¹ | - |
| `m1`, `m2` | `Float64` | 两条传播线的夸克质量 | fm⁻¹ | - |
| `μ1`, `μ2` | `Float64` | 两条传播线的化学势 | fm⁻¹ | - |
| `T` | `Float64` | 温度 | fm⁻¹ | - |
| `Φ`, `Φbar` | `Float64` | Polyakov 圈参数 | 无量纲 | - |
| `ξ` | `Float64` | 各向异性参数 | 无量纲 | - |
| `rtol` | `Float64` | 相对误差容限 | - | `1e-3` |
| `atol` | `Float64` | 绝对误差容限 | - | `0.0` |

#### 返回值

返回元组 `(real_correction, imag_correction)`，分别对应 B₀ 修正项的实部和虚部（单位：fm）。

**重要**：需与各向同性 B₀ 相加：
```julia
B0_iso_real, B0_iso_imag = OneLoopIntegrals.B0(λ, k, m1, μ1, m2, μ2, T; Φ=Φ, Φbar=Φbar)
ΔB0_real, ΔB0_imag = B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ)
B0_total = (B0_iso_real + ΔB0_real, B0_iso_imag + ΔB0_imag)
```

#### 物理意义

B₀ 积分在极化函数中的作用：
```math
\Pi \propto A_1 + A_2 + [k^2 - (k_0 + \mu_1 - \mu_2)^2 + (m_1 \mp m_2)^2] \times B_0
```

各向异性修正通过分布函数导数项传播到 B₀，影响极化函数的虚部（吸收/产生过程）和实部（质量重整化）。

#### 数值实现要点

- k=0 和 k>0 使用不同的积分策略（与 `B0` 相同）
- 虚部计算考虑奇点区间的变化（各向异性影响奇点位置）
- ξ=0 时应返回 (0.0, 0.0)
- 线性性质：`B0_correction(2ξ) / B0_correction(ξ) ≈ 2.0`（小 ξ 时）

#### 使用示例

```julia
using Pkg
Pkg.activate(".")

include("src/relaxtime/OneLoopIntegrals.jl")
include("src/relaxtime/OneLoopIntegralsAniso.jl")

using .OneLoopIntegrals: B0
using .OneLoopIntegralsCorrection: B0_correction

# 物理参数
λ = 0.45
k = 0.30
m1 = 0.24
m2 = 0.38
μ1 = 0.12
μ2 = -0.05
T = 0.17
Φ = 0.15
Φbar = 0.15
ξ = 0.1

# 计算各向同性 B₀
B0_iso = B0(λ, k, m1, μ1, m2, μ2, T; Φ=Φ, Φbar=Φbar)

# 计算一阶修正
ΔB0 = B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ)

# 总结果
B0_total = (B0_iso[1] + ΔB0[1], B0_iso[2] + ΔB0[2])

println("B₀ (各向同性):")
println("  实部: ", B0_iso[1], " fm")
println("  虚部: ", B0_iso[2], " fm")
println("ΔB₀ (一阶修正, ξ=", ξ, "):")
println("  实部: ", ΔB0[1], " fm")
println("  虚部: ", ΔB0[2], " fm")
println("B₀ (总计):")
println("  实部: ", B0_total[1], " fm")
println("  虚部: ", B0_total[2], " fm")
println("修正比例:")
println("  实部: ", 100*ΔB0[1]/B0_iso[1], " %")
println("  虚部: ", 100*ΔB0[2]/B0_iso[2], " %")
```

**典型输出**：
```
B₀ (各向同性):
  实部: -0.05327 fm
  虚部: 0.00123 fm
ΔB₀ (一阶修正, ξ=0.1):
  实部: -0.00181 fm
  虚部: 0.00004 fm
B₀ (总计):
  实部: -0.05508 fm
  虚部: 0.00127 fm
修正比例:
  实部: 3.4 %
  虚部: 3.2 %
```

---

## 便捷工具函数

### `build_default_nodes_weights()`

生成默认的积分节点和权重，适用于标准精度计算。

#### 函数签名

```julia
build_default_nodes_weights() -> (nodes_cosθ, weights_cosθ, 
                                    nodes_p_inf, weights_p_inf,
                                    nodes_p_Λ, weights_p_Λ)
```

#### 返回值

| 变量 | 说明 | 区间 | 节点数 |
|------|------|------|--------|
| `nodes_cosθ` | 角度积分节点 | [-1, 1] | 32 |
| `weights_cosθ` | 角度积分权重 | - | 32 |
| `nodes_p_inf` | 动量节点（扩展范围） | [0, 20 fm⁻¹] | 64 |
| `weights_p_inf` | 动量权重（扩展范围） | - | 64 |
| `nodes_p_Λ` | 动量节点（截断到Λ） | [0, Λ_inv_fm] | 64 |
| `weights_p_Λ` | 动量权重（截断到Λ） | - | 64 |

#### 使用示例

```julia
include("src/integration/GaussLegendre.jl")
using .GaussLegendre: build_default_nodes_weights

nodes_cosθ, weights_cosθ, nodes_p, weights_p, _, _ = build_default_nodes_weights()

# 用于 A_aniso
A_result = A_aniso(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p, nodes_cosθ, weights_cosθ)
```

---

## 性能优化建议

### 1. 节点复用

积分节点生成较耗时，应在循环外预生成：

```julia
# 好的做法：节点复用
nodes_p, weights_p = gauleg(0.0, 20.0, 64)
nodes_cosθ, weights_cosθ = gauleg(-1.0, 1.0, 32)

for ξ in [0.0, 0.1, 0.2, 0.3]
    A_result = A_aniso(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p, nodes_cosθ, weights_cosθ)
    println("ξ = ", ξ, ", A = ", A_result)
end

# 坏的做法：重复生成节点
for ξ in [0.0, 0.1, 0.2, 0.3]
    nodes_p, weights_p = gauleg(0.0, 20.0, 64)  # 每次循环都重新生成！
    A_result = A_aniso(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p, nodes_cosθ, weights_cosθ)
end
```

### 2. 各向同性项复用

扫描多个 ξ 值时，各向同性项 `A(m, μ, T, Φ, Φbar, nodes_p, weights_p)` 只需计算一次：

```julia
# 优化：A_iso 只计算一次
A_iso = A(m, μ, T, Φ, Φbar, nodes_p, weights_p)

results = []
for ξ in [0.05, 0.10, 0.15, 0.20]
    ΔA = A_correction(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p)
    push!(results, (ξ, A_iso + ΔA))
end
```

### 3. 积分精度权衡

- **快速计算**：`nodes_p` 32个点，`nodes_cosθ` 16个点（相对误差~1%）
- **标准计算**：`nodes_p` 64个点，`nodes_cosθ` 32个点（相对误差~0.1%）
- **高精度计算**：`nodes_p` 128个点，`nodes_cosθ` 64个点（相对误差~0.01%）

### 4. 适用范围选择

| ξ 范围 | 推荐方法 | 理由 |
|--------|---------|------|
| 0-0.2 | `A + A_correction` | 一阶近似足够，性能最优 |
| 0.2-0.4 | 两者均可 | 误差<5%，根据需求选择 |
| 0.4-1.0 | `A_aniso` | 一阶近似失效，需完整计算 |

---

## 验证与测试

### 边界条件检验

```julia
# 1. ξ=0 时，修正项应为零
ΔA = A_correction(m, μ, T, Φ, Φbar, 0.0, nodes_p, weights_p)
@assert abs(ΔA) < 1e-10  # 应几乎为零

# 2. ξ=0 时，A_aniso 应等于 A_iso
A_iso = A(m, μ, T, Φ, Φbar, nodes_p, weights_p)
A_aniso_zero = A_aniso(m, μ, T, Φ, Φbar, 0.0, nodes_p, weights_p, nodes_cosθ, weights_cosθ)
@assert abs(A_aniso_zero - A_iso) / A_iso < 1e-3  # 相对误差<0.1%
```

### 线性响应检验

```julia
# 小 ξ 时，修正项应线性于 ξ
ξ1 = 0.05
ξ2 = 0.10
ΔA1 = A_correction(m, μ, T, Φ, Φbar, ξ1, nodes_p, weights_p)
ΔA2 = A_correction(m, μ, T, Φ, Φbar, ξ2, nodes_p, weights_p)
ratio = ΔA2 / ΔA1
@assert abs(ratio - 2.0) < 0.05  # 比值应接近 2.0（5%误差容忍）
```

### 收敛性测试

```julia
# 增加节点数，结果应收敛
for n_points in [32, 64, 128, 256]
    nodes, weights = gauleg(0.0, 20.0, n_points)
    A_val = A_correction(m, μ, T, Φ, Φbar, ξ, nodes, weights)
    println("n=", n_points, " -> A_correction=", A_val)
end
# 预期：n≥64 后结果变化<1%
```

---

## 常见问题

### Q1: `A + A_correction` 和 `A_aniso` 结果不一致？

**A**: 正常现象。一阶近似只保留 ξ 的线性项，而 `A_aniso` 包含所有高阶项。差异大小取决于 ξ：
- ξ=0.1：差异通常<1%
- ξ=0.3：差异约2-5%
- ξ>0.5：差异可达10%以上

### Q2: 为什么 `A_correction` 可能为负值？

**A**: 各向异性沿 z 轴时，垂直方向的分布被"压缩"，某些动量区域的有效态密度降低，导致负修正。总的 `A_total` 仍为正。

### Q3: `B0_correction` 的虚部什么时候非零？

**A**: 当 k² 和 λ 满足特定关系，使得积分奇点落在物理区域内时。各向异性会轻微改变奇点位置，但通常虚部修正很小（<5%）。

### Q4: 如何选择动量积分上限？

**A**: 默认 20 fm⁻¹ 对大多数情况足够。如果：
- T 很高（>0.4 fm⁻¹）：增加到 30 fm⁻¹
- μ 很大（>0.5 fm⁻¹）：增加到 25 fm⁻¹
- 检查收敛性：逐步增大上限，看结果是否稳定

---

## 相关文档

- **公式推导**：`docs/reference/formula/relaxtime/integrals/OneLoopIntegral_A_Aniso_Collinear.md`
- **各向同性基础**：`docs/api/integrals/OneLoopIntegrals.md`
- **分布函数**：`docs/reference/formula/pnjl/distribution/PNJL_夸克有效分布函数_动量各向异性.md`
- **测试总结**：`test/test_oneloopintegrals_aniso_summary.md`
- **B0修正详细分析**：`test/test_b0_correction_summary.md`

---

## 更新日志

- **2025-12-24**: 增强 `tilde_B0_correction_k_positive` 函数，添加聚簇 GL 积分策略，显著提升 A±B 根处的数值稳定性
- **2025-11-14**: 初始版本，添加 `A_correction`、`A_aniso`、`B0_correction` 的完整 API 文档

---

## 数值积分策略（v2.0 新增）

### 背景

在计算 `tilde_B0_correction_k_positive` 时，被积函数在 A±B=0 的根处存在对数奇异性，导致标准 QuadGK 积分器在某些参数组合下精度下降。为解决此问题，我们实现了多种数值积分策略。

### 可用策略

| 策略 | 枚举值 | 说明 | 推荐场景 |
|------|--------|------|----------|
| QuadGK | `STRATEGY_QUADGK` | 原始自适应积分 | 快速估算、调试 |
| 区间分割 GL | `STRATEGY_INTERVAL_GL` | 在根处分割区间 + 标准 GL | 中等精度需求 |
| 聚簇 GL | `STRATEGY_CLUSTER_GL` | 在根处分割区间 + tanh 映射聚簇 GL | 高精度需求 |
| 混合策略 | `STRATEGY_HYBRID` | 根据奇点位置自适应选择变换 | **默认**，最优性能 |

### 混合策略原理

STRATEGY_HYBRID 根据奇点在区间中的位置选择最优变换：

1. **单侧奇点（左端）**：使用幂次变换 `t = u^(1/α)` 聚簇于左端
2. **单侧奇点（右端）**：使用幂次变换 `t = 1 - (1-u)^(1/α)` 聚簇于右端
3. **双侧奇点**：使用 DE (tanh-sinh) 变换，两端都有聚集
4. **无奇点**：使用 `power_left` 聚簇于低能端

无根情况下使用 `power_left` 的原因：被积函数包含 `p = sqrt(E² - m²)`，在 E→m（低能端）时变化最剧烈，所有测试参数组合都验证了这一点。

### 聚簇 GL 原理

聚簇 GL 使用 tanh 映射将标准 GL 节点聚集到区间端点附近：

```math
x = \frac{a+b}{2} + \frac{b-a}{2} \cdot \frac{\tanh(\beta u)}{\tanh(\beta)}, \quad u \in [-1, 1]
```

其中 β 控制聚簇程度：
- β=1: 轻微聚簇
- β=4: 中等聚簇
- β=8: 强聚簇（推荐）
- β>16: 过度聚簇，可能导致精度下降

### 使用示例

```julia
using .OneLoopIntegralsCorrection: tilde_B0_correction_k_positive,
    STRATEGY_QUADGK, STRATEGY_INTERVAL_GL, STRATEGY_CLUSTER_GL

# 默认调用（使用 STRATEGY_CLUSTER_GL）
result = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ)

# 指定策略
result_quadgk = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=STRATEGY_QUADGK)

# 自定义聚簇参数
result_custom = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=STRATEGY_CLUSTER_GL,
    cluster_beta=8.0,  # 聚簇强度
    cluster_n=64       # 每区间节点数
)

# 获取诊断信息
real_part, imag_part, diag = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    diagnostics=true)

println("策略: ", diag.strategy)
println("找到的根数: ", diag.n_roots)
println("根位置: ", diag.roots)
println("区间数: ", diag.n_intervals)
println("耗时: ", diag.elapsed_ms, " ms")
```

### 精度对比

典型参数下的相对误差（使用默认 n=32）：

| 参数组合 | CLUSTER_GL | HYBRID |
|----------|------------|--------|
| λ=-1, k=0.01 (有根) | 9.8e-4 | 3.8e-6 |
| λ=-1, k=0.1 (有根) | 5.8e-3 | 9.6e-6 |
| λ=-0.5, k=0.05 (无根) | 4.5e-2 | 2.9e-14 |
| λ=0.5, k=0.1 (无根) | 6.8e-2 | 1.6e-13 |

**结论**：HYBRID 策略在所有情况下都能以 n=32 达到 1e-5 精度，是默认推荐选项。

### 默认参数

```julia
DEFAULT_STRATEGY = STRATEGY_HYBRID  # 混合策略（最优性能）
DEFAULT_CLUSTER_BETA = 8.0          # tanh 聚簇参数
DEFAULT_CLUSTER_N = 32              # 每区间节点数（实际使用自适应）
DEFAULT_POWER_ALPHA = 0.35          # 幂次变换参数
DEFAULT_DE_H = 0.15                 # DE 变换步长
```

### 自适应节点数分配

HYBRID 策略根据区间类型自适应分配节点数：

| 情况 | 区间类型 | 节点数 | 说明 |
|------|----------|--------|------|
| 无根 | 单区间 (SING_NONE) | 16 | 使用 power_left 聚簇于低能端 |
| 有根 | 第一个区间 (SING_RIGHT) | 32 | 左端 E=m 处有 p→0 奇异行为 |
| 有根 | 中间区间 (SING_BOTH) | 32 | 两端都有奇点 |
| 有根 | 最后一个区间 (SING_LEFT) | 16 | 只有左端奇点 |

**典型节点总数**：
- 无根情况：16 节点
- 1 根情况：32 + 16 = 48 节点
- 2 根情况：32 + 32 + 16 = 80 节点

### 性能对比

| 策略 | 平均用时 | 精度达标率 | 最大误差 |
|------|----------|------------|----------|
| QuadGK (rtol=1e-3) | 4.3 μs | 低 | 2.9e-2 |
| **HYBRID (自适应)** | **1.1 μs (无根) / 3.7 μs (有根)** | **高** | **1.6e-5** |

**说明**：
- 无根情况：HYBRID 使用 16 节点，约 1.1 μs，比 QuadGK 快 4 倍
- 有根情况：HYBRID 使用 48 节点，约 3.7 μs，精度远高于 QuadGK
- 对于典型物理参数（m ≥ 0.3 GeV），精度达标率 100%
- 预计算的标准 GL 节点和变换系数避免了重复计算

### 根的解析求解

A±B=0 的根可以通过二次方程解析求解：

```
(k² - λ²)E² - λCE - (k²m² + C²/4) = 0
```

其中 `C = λ² + m² - m'² - k²`。这比数值扫描快约 18 倍。

### 诊断输出结构

```julia
struct IntegrationDiagnostics
    strategy::IntegrationStrategy  # 使用的策略
    n_roots::Int                   # 找到的 A±B=0 根数
    roots::Vector{Float64}         # 根位置
    n_intervals::Int               # 分割的区间数
    intervals::Vector{Tuple{Float64,Float64}}  # 区间列表
    real_part::Float64             # 计算的实部
    imag_part::Float64             # 计算的虚部
    elapsed_ms::Float64            # 耗时（毫秒）
end
```
