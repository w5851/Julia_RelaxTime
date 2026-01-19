---
title: 介子质量复数 p0 接口实现
archived: true
original: docs/dev/active/介子质量复数p0接口草案.md
archived_date: 2026-01-19
---

以下为原始内容（保留，以便审阅与历史参考）：

***

# 介子质量复数 $p_0$ 接口草案（函数签名 + 伪代码 + 参数说明）

> 目标：在不改变现有实数积分路径的前提下，支持 $p_0 = M + i\Gamma/2$ 的极点求解。
> 关键原则：**宽度 $\Gamma$（或 `gamma`）不进入 B0 的积分路径**，而是在极化函数组合式中与 $B_{\text{re}}, B_{\text{im}}$ 发生耦合。

## 1. 背景依据（与 Fortran 版本一致）
Fortran 中的实现方式是：
- 以实数 `k0` 计算 B0 的实部/虚部；
- 通过显式组合公式把 `ga`（宽度）耦合进 $\Pi_{\text{re}}, \Pi_{\text{im}}$；
- **没有对复数能量直接积分**。

因此 Julia 侧若扩展复数 $p_0$，建议采用**“拆分 $p_0$”**的方案：
- $k0 = \Re(p_0)$
- $\gamma = 2\Im(p_0)$
- 依照下述组合公式计算 $\Pi$。

---

## 2. 函数签名（草案）

### 2.1 现有实数极化函数（保持不变）
```
polarization_aniso(channel::Symbol,
                   k0::Float64, k_norm::Float64,
                   m1::Float64, m2::Float64,
                   μ1::Float64, μ2::Float64,
                   T::Float64, Φ::Float64, Φbar::Float64,
                   ξ::Float64, A1::Float64, A2::Float64,
                   num_s::Int)
    -> (Π_real::Float64, Π_imag::Float64)
```

### 2.2 新增：带宽度的极化函数
```
polarization_with_width(channel::Symbol,
                        k0::Float64, gamma::Float64, k_norm::Float64,
                        m1::Float64, m2::Float64,
                        μ1::Float64, μ2::Float64,
                        T::Float64, Φ::Float64, Φbar::Float64,
                        ξ::Float64, A1::Float64, A2::Float64,
                        num_s::Int)
    -> (Π_real::Float64, Π_imag::Float64)
```

### 2.3 可选：复数能量包装接口
```
polarization_complex(channel::Symbol,
                     p0::ComplexF64, k_norm::Float64,
                     m1::Float64, m2::Float64,
                     μ1::Float64, μ2::Float64,
                     T::Float64, Φ::Float64, Φbar::Float64,
                     ξ::Float64, A1::Float64, A2::Float64,
                     num_s::Int)
    -> (Π_real::Float64, Π_imag::Float64)
```

---

## 3. 伪代码（核心流程）

### 3.1 `polarization_with_width`
```
function polarization_with_width(channel, k0, gamma, k_norm,
                                 m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, num_s):
    # 1) 实数路径：计算 B0 的实部/虚部
    (B_re, B_im) = compute_B0_parts(k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, num_s)

    # 2) 计算 λ 与组合项
    λ = k0 + μ1 - μ2

    # channel-dependent mass term
    if channel == :P:
        mass_term = (m1 - m2)^2
    else if channel == :S:
        mass_term = (m1 + m2)^2
    else:
        throw error

    C = mass_term - λ^2 + (gamma^2)/4

    # 3) 组合极化函数实部/虚部
    Π_real = factor * (A1 + A2 + C * B_re - gamma * λ * B_im)
    Π_imag = factor * (C * B_im + gamma * λ * B_re)

    return (Π_real, Π_imag)
```

> 说明：`factor` 与现有 `polarization_aniso` 中的系数保持一致。

### 3.2 `polarization_complex`
```
function polarization_complex(channel, p0, k_norm, ...):
    k0 = real(p0)
    gamma = 2 * imag(p0)
    return polarization_with_width(channel, k0, gamma, k_norm, ...)
```

---

## 4. 参数说明

### 通用参数
- `channel::Symbol`：
  - `:P` 赝标量通道（π、K、η、η′ 等）
  - `:S` 标量通道（σ、σ′ 等）

- `k0::Float64`：实数能量分量（$\Re(p_0)$）。
- `gamma::Float64`：衰变宽度（$\Gamma$），对应 $p_0 = k0 + i\Gamma/2$。
- `k_norm::Float64`：三动量模长 $|\vec{k}|$。
- `m1, m2::Float64`：两味夸克动力学质量。
- `μ1, μ2::Float64`：两味夸克化学势。
- `T::Float64`：温度。
- `Φ, Φbar::Float64`：Polyakov loop。
- `ξ::Float64`：各向异性参数。
- `A1, A2::Float64`：A 函数预计算值。
- `num_s::Int`：奇异夸克数量（0 或 1）。

### 返回值
- `(Π_real, Π_imag)`：极化函数实部与虚部，供传播子与极点方程使用。

---

## 5. 设计说明（为何不直接用复数积分）
- 复数 $p_0$ 直接进入积分会触发解析延拓与分支选择问题，数值实现复杂且稳定性差。
- Fortran 版本明确采用“实数积分 + 宽度耦合”的结构，验证过可用。
- 该方案能与现有缓存、B0 计算流程保持一致，只需在上层组合公式中引入 $\gamma$。

---

## 6. 实施注意点（未来落地时）
- 若引入 `gamma`，缓存键需要新增 `gamma` 字段（不建议直接存 `Complex`）。
- 极点求解层应同时求解 `Re` 与 `Im` 方程，或解 `fvec=(1-4KΠ_re, Π_im)` 形式。
- 文档中需明确：`B0` 仍按实数路径计算，`gamma` 仅在组合式出现。