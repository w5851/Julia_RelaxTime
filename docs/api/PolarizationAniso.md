# PolarizationAniso 模块 API 文档

## 模块概述

`PolarizationAniso` 模块提供 PNJL 模型中赝标量 (`:P`) 与标量 (`:S`) 介子极化函数在**动量各向异性**条件下的数值计算。实现依赖于 `OneLoopIntegrals` 模块中的单圈积分 `B0` 与单线积分 `A`，以及 `OneLoopIntegralsAniso` 模块中的各向异性修正 `B0_correction`，并支持复用预计算的 `A` 值以减少重复工作量。

## 依赖

- `Constants_PNJL` —— 提供颜色数等常量
- `OneLoopIntegrals` —— 提供 `B0` 与 `A` 的数值实现
- `OneLoopIntegralsAniso` —— 提供 `B0_correction` 各向异性修正

## 单位约定

同项目其它模块一致，全部使用自然单位制 (ℏ = c = 1)：

- `k0`, `k_norm`, `m1`, `m2`, `μ1`, `μ2`, `T` 采用 fm⁻¹
- `Φ`, `Φbar` 为无量纲 Polyakov 圈变量
- `ξ` 为无量纲各向异性参数，取值范围 `[-1, 1]`，`0` 表示各向同性
- 返回的极化函数实部与虚部同样以 fm 单位计

---

## API 参考

### `polarization_aniso(channel, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1_value, A2_value, num_s_quark)`

计算给定通道 (`:P` 赝标量或 `:S` 标量) 下的各向异性极化函数，返回 `(ReΠ, ImΠ)`。

#### 函数签名

```julia
polarization_aniso(channel::Symbol, k0::Float64, k_norm::Float64,
                   m1::Float64, m2::Float64,
                   μ1::Float64, μ2::Float64,
                   T::Float64, Φ::Float64, Φbar::Float64,
                   ξ::Float64,
                   A1_value::Float64, A2_value::Float64,
                   num_s_quark::Int) -> NTuple{2, Float64}
```

#### 参数

| 参数 | 类型 | 说明 | 单位 |
|------|------|------|------|
| `channel` | `Symbol` | `:P` (pseudoscalar) 或 `:S` (scalar) | — |
| `k0` | `Float64` | 外部能量 | fm⁻¹ |
| `k_norm` | `Float64` | 外部三动量模长 | fm⁻¹ |
| `m1`, `m2` | `Float64` | 两条传播线的质量 | fm⁻¹ |
| `μ1`, `μ2` | `Float64` | 两条传播线的化学势 | fm⁻¹ |
| `T` | `Float64` | 温度 | fm⁻¹ |
| `Φ`, `Φbar` | `Float64` | Polyakov 圈变量 | 无量纲 |
| `ξ` | `Float64` | 各向异性参数，范围 `[-1, 1]`；`ξ = 0` 表示各向同性，`ξ > 0` 表示沿动量方向拉伸，`ξ < 0` 表示垂直方向拉伸 | 无量纲 |
| `A1_value`, `A2_value` | `Float64` | 对应传播线的单线积分 `A` 结果，可由 `OneLoopIntegrals.A` 预先计算 | fm |
| `num_s_quark` | `Int` | strange 夸克数开关；`0`或`2` 表示直接使用 `B0`，`1` 表示对 `k0` 与 `-k0` 的 `B0` 结果取平均以恢复对称性 | — |

#### 返回值

返回 `(real_part, imag_part)` 二元组，对应极化函数的实部与虚部。

#### 数值实现要点

- 预先调用 `OneLoopIntegrals.A` 生成 `A1_value`、`A2_value` 可避免重复积分。
- 计算中依赖 `B0` 与 `B0_correction`，因此与它们一样可通过关键字参数控制积分精度。
- `num_s_quark = 1` 时，内部会对 `B0(λ, k)` 与 `B0(-k0 + μ1 - μ2, k)` 的结果取平均，以缓解奇点导致的对称性破坏；其它取值留作扩展。
- `channel` 仅接受 `:P` 与 `:S`，传入其它符号会抛出 `ArgumentError`。
- 当 `ξ > EPS_SEGMENT` 时（`EPS_SEGMENT` 定义于 `OneLoopIntegrals` 模块），会计算并应用各向异性修正。

#### 各向异性参数 ξ 的物理意义

各向异性参数 `ξ` 描述动量空间分布的各向异性程度：

- **ξ = 0**：各向同性情况，退化为标准的 `polarization` 函数
- **ξ > 0**：沿某一特定方向（通常为外场或碰撞轴方向）的动量分布被拉伸
- **ξ < 0**：垂直于特定方向的动量分布被拉伸
- **典型取值范围**：`[-1, 1]`，极端情况下可能超出此范围

#### 使用示例

```julia
include("src/relaxtime/PolarizationAniso.jl")
using .PolarizationAniso
include("src/relaxtime/OneLoopIntegrals.jl")
using .OneLoopIntegrals: A, gauleg

nodes, weights = gauleg(0.0, 20.0, 128)
A1 = A(0.3, 0.1, 0.15, 0.2, 0.2, nodes, weights)
A2 = A(0.35, 0.05, 0.15, 0.2, 0.2, nodes, weights)

# 各向同性情况 (ξ = 0)
ΠP_iso = polarization_aniso(:P, 0.4, 0.2, 0.3, 0.35, 0.1, 0.05, 0.15, 0.2, 0.2, 0.0, A1, A2, 0)

# 各向异性情况 (ξ = 0.5)
ΠP_aniso = polarization_aniso(:P, 0.4, 0.2, 0.3, 0.35, 0.1, 0.05, 0.15, 0.2, 0.2, 0.5, A1, A2, 0)

# 强各向异性 (ξ = -0.8)
ΠS_aniso = polarization_aniso(:S, 0.4, 0.2, 0.3, 0.35, 0.1, 0.05, 0.15, 0.2, 0.2, -0.8, A1, A2, 1)
```

---

## 数值与性能建议

- `A` 的计算仅依赖质量、化学势、温度及 Polyakov 圈，适合在外层缓存以复用。
- 调整 `B0` 与 `B0_correction` 调用时的 `rtol`、`atol` 可在极化函数需要更高精度时收紧误差。
- `num_s_quark = 1` 可用于含奇异夸克的通道以保持 `k0` ↔ `-k0` 对称性；其它情形通常维持 `0` 即可。
- 对多组 (`k0`, `k_norm`, `ξ`) 参数进行扫描时，可预先生成 `A` 并共享。
- 当 `ξ` 接近 `0` 时（小于 `EPS_SEGMENT`），各向异性修正会被自动跳过以节省计算时间。
- 各向异性修正的计算成本较高，建议在需要时才使用非零 `ξ` 值。

## 与 Polarization 模块的关系

- `PolarizationAniso` 是 `Polarization` 模块的扩展版本，增加了各向异性效应。
- 当 `ξ = 0` 且 `ξ <= EPS_SEGMENT` 时，`polarization_aniso` 的结果应与 `polarization` 一致。
- 推荐在各向同性研究中使用 `Polarization` 模块以获得更好的性能，在需要考虑各向异性效应时使用 `PolarizationAniso` 模块。

## 更新日志

- **v0.1.0 (2025-11-10)**：初始撰写 PolarizationAniso API 文档，覆盖函数签名、参数说明、各向异性物理意义与示例。
