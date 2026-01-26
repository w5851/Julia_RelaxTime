# 调查总结：Fortran vs Julia 弛豫时间差异

**日期**: 2026-01-26
**状态**: 已找到关键线索，需要进一步验证

---

## 🎯 核心发现

### 1. 数密度一致 ✅
- Fortran: ρ_u = 1.681 fm⁻³, ρ_s = 1.538 fm⁻³
- Julia: ρ_u = 1.694 fm⁻³, ρ_s = 1.550 fm⁻³
- **差异 < 1%**

### 2. A 值一致 ✅
- Fortran: A_u = -4.972 fm⁻², A_s = -5.166 fm⁻²
- Julia: A_u = -4.964 fm⁻², A_s = -5.158 fm⁻²
- **差异 < 0.2%**

### 3. G 值有差异 ❌
- Fortran: arrG(u) = 0.000714 fm², arrG(s) = 0.01891 fm²
- Julia: G_u = 0.015282 fm⁻³, G_s = 0.4046 fm⁻³
- **需要理解单位和定义**

---

## 🔍 关键观察（用户提供）

根据文档 `EffectiveCoupling_K_FromA.md`：

```
K₀± = G ∓ (1/3)K(2G^μ + G^s)
```

这里的 `K·G^μ` 是一个整体项。

**推论**：
- Fortran 计算的 `arrG` 应该是 `K·G^f`（整体，单位 fm²）
- Julia 计算的 `G_f` 应该是纯粹的 `G^f`（无量纲或 fm⁻³）
- **关系**：`arrG (Fortran) = K_f × G^f (Julia)`

---

## 📊 数值验证

### 预期关系

```
arrG(u) = K_f × G^f
```

其中：
- `K_f = 0.00519 fm⁵` (Fortran)
- `arrG(u) = 0.000714 fm²` (Fortran)

因此：
```
G^f = arrG(u) / K_f = 0.000714 / 0.00519 = 0.1376 fm⁻³
```

### Julia 的计算

```
G^f (Julia) = -3 / (4π²) × m_u × A_u
            = -0.07598 × 0.040510 × (-4.964141)
            = 0.015282 fm⁻³
```

### 差异

```
G^f (Fortran) / G^f (Julia) = 0.1376 / 0.015282 = 9.0
```

**差了 9 倍！**

---

## 🎓 可能的原因

### 1. A 的归一化不同

**Fortran**：
```fortran
arrA(:) = arrA(:)*4d0  ! 最后乘以 4
```

**Julia**：
```julia
A_f = 16π² × [积分]  ! 可能包含 16π² 因子
```

**如果**：
- Fortran 的 A 定义中缺少某个因子
- 或者 Julia 的 A 定义中多了某个因子
- 这会导致 G^f 的差异

### 2. G^f 的定义不同

**标准定义**（从文档）：
```
G^f = -N_c/(4π²) · m_f · A_f
```

**Fortran 实现**：
```fortran
arrG = - m * A * (3 * K_f) / (4π²)
     = [- 3/(4π²) · m · A] × K_f
     = G^f × K_f
```

**Julia 实现**：
```julia
G_f = -3 / (4π²) * m_f * A_f
```

**公式相同，但如果 A 的定义不同，结果就会不同。**

### 3. 单位约定不同

**可能性**：
- Fortran 使用自然单位制，某些因子被吸收
- Julia 使用显式单位，所有因子都明确
- 需要检查两者的单位约定

---

## 📝 下一步行动

### 立即需要做的

1. **检查 A 函数的完整定义**：
   - Fortran: `quantity.f90` 中的 A 计算
   - Julia: `OneLoopIntegrals.jl` 中的 A 函数
   - 对比归一化因子

2. **检查 G^f 的使用方式**：
   - Fortran: arrG 如何在传播子中使用
   - Julia: G_f 如何在 `calculate_effective_couplings` 中使用
   - 验证单位是否匹配

3. **检查 K 值的计算**：
   - 从 arrG 和 G_f 计算出的 K 值
   - 对比是否一致

### 验证策略

1. **假设 Fortran 正确**：
   - 修改 Julia 的 `calculate_G_from_A`
   - 使其返回 `K·G^f` 而不是 `G^f`
   - 重新计算弛豫时间

2. **假设 Julia 正确**：
   - 修改 Fortran 的 `arrG` 计算
   - 移除 `K_f` 因子
   - 重新编译运行

3. **对比中间结果**：
   - 计算 K 值（有效耦合系数）
   - 计算 M²（散射振幅）
   - 计算 σ(s)（截面）
   - 找出第一个出现差异的地方

---

## 🎯 预期结果

如果我们的分析正确：

1. **修正后**，Fortran 和 Julia 的 G 值应该一致
2. **修正后**，K 值应该一致
3. **修正后**，弛豫时间应该接近（差异 < 10%）

**关键**：需要确认 `arrG` 和 `G_f` 的正确定义和关系。

---

## 📚 相关文档

- `EffectiveCoupling_K_FromA.md`: 有效耦合系数的公式
- `FORTRAN_JULIA_A_G_COMPARISON.md`: A 和 G 值的详细对比
- `ROOT_CAUSE_FOUND.md`: 初步分析（部分结论需要修正）
- `UNIT_ANALYSIS.md`: 单位分析
- `FINAL_VERIFICATION.md`: 数值验证

---

*文档时间: 2026-01-26*
*状态: 已找到 G 值差异（9倍），需要检查 A 的归一化和 G^f 的定义*
