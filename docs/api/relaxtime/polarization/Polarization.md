# Polarization 模块 API 文档

## 模块概述

`Polarization` 模块提供 PNJL 模型中赝标量 (`:P`) 与标量 (`:S`) 介子极化函数的数值计算。实现依赖于 `OneLoopIntegrals` 模块中的单圈积分 `B0` 与单线积分 `A`，并支持复用预计算的 `A` 值以减少重复工作量。

## 依赖

- `Constants_PNJL` —— 提供颜色数等常量
- `OneLoopIntegrals` —— 提供 `B0` 与 `A` 的数值实现

## 单位约定

同项目其它模块一致，全部使用自然单位制 (ℏ = c = 1)：

- `k0`, `k_norm`, `m1`, `m2`, `μ1`, `μ2`, `T` 采用 fm⁻¹
- `Φ`, `Φbar` 为无量纲 Polyakov 圈变量
- 返回的极化函数实部与虚部同样以 fm 单位计

---

## API 参考

### `polarization(channel, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, A1_value, A2_value, num_s_quark)`

计算给定通道 (`:P` 赝标量或 `:S` 标量) 下的极化函数，返回 `(ReΠ, ImΠ)`。

#### 函数签名

```julia
polarization(channel::Symbol, k0::Float64, k_norm::Float64,
             m1::Float64, m2::Float64,
             μ1::Float64, μ2::Float64,
             T::Float64, Φ::Float64, Φbar::Float64,
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
| `A1_value`, `A2_value` | `Float64` | 对应传播线的单线积分 `A` 结果，可由 `OneLoopIntegrals.A` 预先计算 | fm |
| `num_s_quark` | `Int` | strange 夸克数开关；`0`或`2` 表示直接使用 `B0`，`1` 表示对 `k0` 与 `-k0` 的 `B0` 结果取平均以恢复对称性 | — |

#### 返回值

返回 `(real_part, imag_part)` 二元组，对应极化函数的实部与虚部。

#### 数值实现要点

- 预先调用 `OneLoopIntegrals.A` 生成 `A1_value`、`A2_value` 可避免重复积分。
- 计算中依赖 `B0`，因此与 `B0` 一样可通过其关键字参数控制积分精度。
- `num_s_quark = 1` 时，内部会对 `B0(λ, k)` 与 `B0(-k0 + μ1 - μ2, k)` 的结果取平均，以缓解奇点导致的对称性破坏；其它取值留作扩展。
- `channel` 仅接受 `:P` 与 `:S`，传入其它符号会抛出 `ArgumentError`。

#### 使用示例

```julia
include("src/relaxtime/Polarization.jl")
using .Polarization
include("src/relaxtime/OneLoopIntegrals.jl")
using .OneLoopIntegrals: A, gauleg

nodes, weights = gauleg(0.0, 20.0, 128)
A1 = A(0.3, 0.1, 0.15, 0.2, 0.2, nodes, weights)
A2 = A(0.35, 0.05, 0.15, 0.2, 0.2, nodes, weights)

ΠP = polarization(:P, 0.4, 0.2, 0.3, 0.35, 0.1, 0.05, 0.15, 0.2, 0.2, A1, A2, 0)
ΠS = polarization(:S, 0.4, 0.2, 0.3, 0.35, 0.1, 0.05, 0.15, 0.2, 0.2, A1, A2, 1)
```

---

## 数值与性能建议

- `A` 的计算仅依赖质量、化学势、温度及 Polyakov 圈，适合在外层缓存以复用。
- 调整 `B0` 调用时的 `rtol`、`atol` 可在极化函数需要更高精度时收紧误差。
- `num_s_quark = 1` 可用于含奇异夸克的通道以保持 `k0` ↔ `-k0` 对称性；其它情形通常维持 `0` 即可。
- 对多组 (`k0`, `k_norm`) 参数进行扫频时，可预先生成 `A` 并共享。

## 更新日志

- **v0.2.0 (2025-11-04)**：新增 `num_s_quark` 参数说明及使用建议。
- **v0.1.0 (2025-11-04)**：初始撰写 Polarization API 文档，覆盖函数签名、参数说明与示例。