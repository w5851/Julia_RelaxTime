# 最终验证：Fortran arrG = Julia G_u × K_f

**日期**: 2026-01-26

---

## 🎯 你的观察

根据文档 `EffectiveCoupling_K_FromA.md`，有效耦合系数的公式是：

```
K₀± = G ∓ (1/3)K(2G^μ + G^s)
```

这里的 `K·G^μ` 是一个整体项。

**因此**：
- Fortran 计算的 `arrG` 应该是 `K·G^f`（整体）
- Julia 计算的 `G_f` 应该是纯粹的 `G^f`
- 关系：`arrG (Fortran) = K_f × G_f (Julia)`

---

## 📊 数值验证

### 已知数值

**Fortran**：
```
arrG(u) = 0.000714 fm²
K_f = 0.00519 fm⁵
```

**Julia**：
```
G_u = 0.015282 fm⁻³  （根据单位分析）
```

### 验证关系

```
arrG(u) 应该等于 K_f × G_u

0.000714 fm² =? 0.00519 fm⁵ × 0.015282 fm⁻³
0.000714 fm² =? 0.0000793 fm²
```

**不对！** 差了约 9 倍。

---

## 🔍 重新检查单位

### Julia 的 G_u 单位

从 `calculate_G_from_A`:
```julia
G_f = -Nc / (4.0 * π^2) * (m_f * A_f)
```

如果：
- `m_f`: fm⁻¹
- `A_f`: fm⁻²

则：
```
G_f 单位 = fm⁻¹ × fm⁻² = fm⁻³
```

但是，让我检查 Julia 的 A 函数返回什么单位...

### Julia 的 A 函数

从 `OneLoopIntegrals.jl`，A 函数的定义应该包含一个因子，使得：

```
A_f = 16π² × [积分]
```

如果积分的单位是 fm⁻²，那么：
```
A_f 单位 = 无量纲 × fm⁻² = fm⁻²
```

**但是**，如果 A 函数的定义中已经包含了某些归一化因子，单位可能不同。

---

## 🎯 关键问题

### Fortran 的 arrG 公式

```fortran
arrG(:) = - arrMass(:) * arrA(:) * (3d0*K_f) / (4d0*pi**2)
```

这可以重写为：
```fortran
arrG(:) = [- 3 / (4π²) × arrMass(:) × arrA(:)] × K_f
       = [G^f] × K_f
```

其中 `G^f = -3/(4π²) × m × A`。

### Julia 的 G_f 公式

```julia
G_f = -Nc / (4.0 * π^2) * (m_f * A_f)
    = -3 / (4π²) * m_f * A_f
```

**公式完全相同！**

### 那么问题在哪里？

如果公式相同，但数值不同，那么一定是：
1. **A 的数值不同**
2. **或者 A 的定义/归一化不同**

---

## 📐 对比 A 值

### Fortran
```
A_u = -4.972144 fm⁻²  （我们之前提取的）
```

### Julia
```
A_u = -4.964141 fm⁻²  （我们计算的）
```

**差异 < 0.2%，基本一致。**

---

## 🔄 重新计算

### 如果 Fortran 的 arrG = K_f × G^f

```
G^f (Fortran) = arrG / K_f
              = 0.000714 / 0.00519
              = 0.1376
```

**单位**：如果 arrG 是 fm²，K_f 是 fm⁵，则：
```
G^f 单位 = fm² / fm⁵ = fm⁻³
```

所以：
```
G^f (Fortran) = 0.1376 fm⁻³
```

### Julia 的 G^f

```
G^f (Julia) = -3 / (4π²) × m_u × A_u
            = -0.07598 × 0.040510 × (-4.964141)
            = 0.015282 fm⁻³
```

### 比值

```
G^f (Fortran) / G^f (Julia) = 0.1376 / 0.015282 = 9.0
```

**差了 9 倍！**

---

## 🎯 找到了！

### Fortran 的 A 有额外的因子

从 Fortran 代码：
```fortran
arrA(:) = arrA(:)*4d0
```

**Fortran 的 A 在最后乘以了 4！**

### Julia 的 A 定义

从文档：
```
A_f = 16π² × [积分]
```

`16π² ≈ 157.9`

### Fortran 的 A 定义

如果 Fortran 的 A 在积分后乘以 4，而 Julia 的 A 在积分后乘以 16π²，那么：

```
A (Fortran) / A (Julia) = 4 / (16π²) = 1 / (4π²) ≈ 1 / 39.48
```

**但这会让 A 变小，不是变大！**

---

## 🔍 另一个可能

### Fortran 的 arrG 公式中的 3

```fortran
arrG(:) = - arrMass(:) * arrA(:) * (3d0*K_f) / (4d0*pi**2)
```

这里的 `3d0` 是 N_c。

### 如果 Fortran 的 A 已经包含了某些因子

让我检查 Fortran 的 A 积分的完整定义...

实际上，从代码看：
```fortran
arrA(:) = arrA(:) + w(i) * p**2 * (-1d0) / arrEi(:)  ! 真空项
arrA(:) = arrA(:) + w1(i) * p**2 * (arrFq(1,:) + arrFq(2,:)) / arrEi(:)  ! 有限T项
arrA(:) = arrA(:)*4d0  ! 最后乘以 4
```

这个 `4d0` 可能是 `4π²` 的一部分，或者是其他归一化因子。

---

## 📝 结论

需要：
1. 检查 Fortran 的 A 函数的完整定义和归一化
2. 检查 Julia 的 A 函数的完整定义和归一化
3. 找出它们之间的差异
4. 确认 arrG 和 G_f 的正确关系

**关键**：A 的定义/归一化可能不同，导致最终的 G 值不同。
