# OneLoopIntegrals 模块 API 文档

## 模块概述

`OneLoopIntegrals` 模块提供有限温度/密度下的单圈极化积分实现，当前包含双传播子积分 `B0` 与单传播子积分 `A` 的数值计算。实现依据 `doc/formula/B0.md` 与 `doc/formula/A.md` 中的推导，支持 PNJL 模型中的有效分布函数，并自动处理三动量为零和非零的两种积分形态。

## 依赖

- `QuadGK.jl` —— 提供自适应一维积分器
- `GaussLegendre` —— 提供高斯-勒让德节点生成（内部工具）
- `PNJLQuarkDistributions` —— 提供夸克/反夸克有效分布函数及其积分
- `Constants_PNJL` —— 提供能量截断 `Λ_inv_fm` 等常量

## 单位约定

模块遵循项目统一的自然单位制 (ℏ = c = 1)：

- 外部参数 `λ`, `k`, `m1`, `m2`, `μ1`, `μ2`, `T` 均以 fm⁻¹ 表示
- 分布函数为无量纲
- 返回值 `(Re, Im)` 分别对应 B₀ 的实部和虚部，以 fm 为单位

---

## API 参考

### `A(m, μ, T, Φ, Φbar, nodes_p, weights_p)`

计算单传播子积分 A。常数项按解析公式截断到 `Λ_inv_fm`，分布函数部分通过高斯-勒让德积分节点进行数值求积。

#### 函数签名

```julia
A(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64,
   nodes_p::Vector{Float64}, weights_p::Vector{Float64}) -> Float64
```

#### 参数

| 参数 | 类型 | 说明 | 单位 |
|------|------|------|------|
| `m` | `Float64` | 传播子质量 | fm⁻¹ |
| `μ` | `Float64` | 化学势 | fm⁻¹ |
| `T` | `Float64` | 温度 | fm⁻¹ |
| `Φ` | `Float64` | Polyakov 圈参数 Φ | 无量纲 |
| `Φbar` | `Float64` | 共轭 Polyakov 圈参数 Φ̄ | 无量纲 |
| `nodes_p` | `Vector{Float64}` | 动量积分节点（通常由 `gauleg` 生成，范围建议覆盖 `[0, p_{max}]`，`p_{max} ≥ Λ_inv_fm`） | fm⁻¹ |
| `weights_p` | `Vector{Float64}` | 对应节点权重 | fm⁻¹ |

#### 返回值

返回 `Float64`，即 A 积分的数值结果。

#### 数值实现要点

- 常数项积分由 `const_integral_term_A` 解析计算，自动使用 `Λ_inv_fm` 作为截断。
- 分布函数项通过节点求积；节点的上限应足够大以覆盖自由粒子尾部（默认推荐 `p_{max} ≈ 15–20 fm⁻¹`）。
- 可通过增大节点数来检验收敛性，不需要额外调整 `rtol`/`atol`。
- `const_integral_term_A` 为内部辅助函数，不在模块导出列表中；如需单独调用，请使用限定名称 `OneLoopIntegrals.const_integral_term_A`。

#### 使用示例

```julia
using RelaxTime
include("src/relaxtime/OneLoopIntegrals.jl")
using .OneLoopIntegrals: A, gauleg

nodes, weights = gauleg(0.0, 20.0, 128)
result = A(0.3, 0.1, 0.15, 0.2, 0.2, nodes, weights)
```

---

### `B0(λ, k, m1, μ1, m2, μ2, T; Φ=0.0, Φbar=0.0, rtol=1.0e-3, atol=0.0)`

计算有限温度/密度下的双传播子单圈积分 B₀，自动组合四个 `\tilde{B}_0` 分量，并根据 `k` 是否接近零选择不同的数值策略。

#### 函数签名

```julia
B0(λ::Float64, k::Float64, m1::Float64, μ1::Float64,
   m2::Float64, μ2::Float64, T::Float64;
   Φ::Float64=0.0, Φbar::Float64=0.0,
   rtol::Float64=1.0e-3, atol::Float64=0.0) -> NTuple{2, Float64}
```

#### 参数

| 参数 | 类型 | 说明 | 单位 | 默认值 |
|------|------|------|------|--------|
| `λ` | `Float64` | Matsubara 频率与化学势的组合量 | fm⁻¹ | — |
| `k` | `Float64` | 外部三动量模长，`k ≈ 0` 时走主值积分分支 | fm⁻¹ | — |
| `m1` | `Float64` | 第一条费米子传播线的质量 | fm⁻¹ | — |
| `μ1` | `Float64` | 第一条传播线的化学势 | fm⁻¹ | — |
| `m2` | `Float64` | 第二条费米子传播线的质量 | fm⁻¹ | — |
| `μ2` | `Float64` | 第二条传播线的化学势 | fm⁻¹ | — |
| `T` | `Float64` | 温度 | fm⁻¹ | — |
| `Φ` | `Float64` | Polyakov 圈参数 Φ | 无量纲 | `0.0` |
| `Φbar` | `Float64` | 共轭 Polyakov 圈参数 Φ̄ | 无量纲 | `0.0` |
| `rtol` | `Float64` | 数值积分相对误差控制 | — | `1e-3` |
| `atol` | `Float64` | 数值积分绝对误差控制 | — | `0.0` |

#### 返回值

返回一个二元组 `(Re, Im)`：

- `Re::Float64` —— B₀ 的实部
- `Im::Float64` —— B₀ 的虚部，虚部仅在满足判据时非零

#### 数值实现要点

- 当 `|k| < 1e-9` (`EPS_K`) 时，函数进入主值积分分支，奇点通过解析处理获得虚部贡献。
- 当 `|k| ≥ 1e-9` 时，函数使用对数形式的解析表达式，并在需要时对奇点区间执行一维积分以获得虚部。
- 截断能量 `Λ_E = √(m² + Λ_inv_fm²)` 自动由模块内部计算。
- 所有积分均采用 `QuadGK.quadgk` 实现，自适应控制误差。

#### 使用示例

```julia
using RelaxTime
include("src/relaxtime/OneLoopIntegrals.jl")
using .OneLoopIntegrals: B0

λ = 0.45
k = 0.30
m1, μ1 = 0.25, 0.10
m2, μ2 = 0.38, -0.05
T = 0.18
Φ = 0.15
Φbar = 0.15

real_part, imag_part = B0(λ, k, m1, μ1, m2, μ2, T; Φ=Φ, Φbar=Φbar)
println("Re(B0) = $real_part, Im(B0) = $imag_part")
```

**高精度计算示例**

```julia
baseline = B0(λ, k, m1, μ1, m2, μ2, T)
high_accuracy = B0(λ, k, m1, μ1, m2, μ2, T; rtol=1e-6)
δRe = abs(baseline[1] - high_accuracy[1])
δIm = abs(baseline[2] - high_accuracy[2])
```

**k → 0 极限验证**

```julia
B0_zero = B0(λ, 0.0, m1, μ1, m2, μ2, T)
B0_small = B0(λ, 1e-4, m1, μ1, m2, μ2, T)
@show B0_zero, B0_small
```

---

## 性能建议

- 对需要高精度的计算，首选调整 `rtol`，`atol` 通常保持零即可。
- 对同一组参数重复调用时，可提前评估所需精度，避免过紧的误差阈值导致无谓的积分细分。
- 当需要批量评估不同的 `λ` 或 `k` 时，建议在外层循环中缓存公共的分布函数参数以减少重复计算。
- 对 `A` 的积分，优先通过增加高斯-勒让德节点数或扩大上限 `p_{max}` 来验证收敛。

## 数值稳定性

- 选择 `λ` 时应避免过接近零的值，因为虚部在 `k = 0` 分支中含有 `1/λ` 项；若物理上必须接近零，可通过外层极限或对称性进行处理。
- 对于较大的 `k`，被积函数中的对数项可能出现数值负值；实现中使用绝对值与自适应积分确保稳定。
- 在 PNJL 模型中，建议保持 `0 ≤ Φ, Φbar ≤ 1` 以避免分布函数的极端尖峰。

## 常见问题

- **返回值为何不是复数？** 当前实现直接返回 `(Re, Im)` 二元组，方便调用侧分别处理实部和虚部；若需要复数，可在调用侧使用 `complex(real_part, imag_part)` 构造。
- **如何验证结果收敛？** 建议在少量代表性参数上降低 `rtol`、提高 `quadgk` 精度，并在测试或基准中比较差异。
- **性能如何评估？** 可使用 `@elapsed` 或 `BenchmarkTools.jl` 对 `B0` 进行多次重复调用，项目测试中也内置了简易性能回归。

## 更新日志

- **v0.2.1 (2025-11-03)**：将 `const_integral_term_A` 调整为内部辅助函数并更新 `A` 的使用说明。
- **v0.2.0 (2025-11-03)**：新增单线积分 `A` 及常数项说明，补充性能建议。
- **v0.1.0 (2025-11-03)**：初始撰写 API 文档，覆盖 `B0` 函数的调用约定、数值策略与性能建议。
