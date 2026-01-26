# Fortran vs Julia 中间量对比

**日期**: 2026-01-26
**目标**: 逐步对比中间量，找出3倍差异的根源

---

## 📊 Julia 计算结果

### 参数设置
- **温度**: T = 300 MeV = 1.520319 fm⁻¹
- **化学势**: μ = 2 MeV = 0.010135 fm⁻¹
- **Polyakov loops**: Φ = 0.840408, Φbar = 0.840412
- **有效质量**: m_u = 0.040510 fm⁻¹ (7.994 MeV), m_s = 1.032295 fm⁻¹ (203.700 MeV)

---

## 1. 数密度 ρ

### Julia 结果
```
ρ_u (Julia)   = 1.693581 fm⁻³
ρ_s (Julia)   = 1.549674 fm⁻³
```

### Fortran 结果
```
ρ_u (Fortran) = 1.680901 fm⁻³
ρ_s (Fortran) = 1.537628 fm⁻³
```

### 对比
```
ρ_u 比值 = 1.0075 (差异 0.75%)
ρ_s 比值 = 1.0078 (差异 0.78%)
```

### 结论
✅ **数密度一致**（差异 < 1%）

**说明**：
- 数密度的计算方法基本一致
- 微小差异可能来自积分精度
- **这不是3倍差异的来源**

---

## 2. 有效耦合常数 K

### Julia 结果

**单圈积分 A 值**：
```
A_u = -4.964141 fm⁻²
A_s = -5.157564 fm⁻²
```

**有效耦合 G 值**：
```
G_u = 0.015282 fm²
G_s = 0.404585 fm²
```

**K 矩阵元素**：
```
K0_plus     = 0.1902
K0_minus    = 0.2037
K123_plus   = 0.2064
K123_minus  = 0.1875
K4567_plus  = 0.1973
K4567_minus = 0.1966
K8_plus     = 0.1943
K8_minus    = 0.1996
```

**物理介子耦合常数**：
```
K_pi    = 0.2064 (K123_plus)
K_K     = 0.1973 (K4567_plus)
K_sigma = 0.1902 (K0_plus)
K_kappa = 0.1966 (K4567_minus)
```

### Fortran 结果

**需要提取**：
1. 从 `coefficients_tmu.f90` 中提取 A_u, A_s 的值
2. 从 `coefficients_tmu.f90` 中提取 G_u, G_s 的值
3. 从 `coefficients_tmu.f90` 中提取 K 矩阵元素

**提取方法**：
- 在 `coefficients_tmu.f90` 中添加 print 语句
- 输出 A_u, A_s, G_u, G_s, K_pi, K_K 等值
- 重新编译运行

---

## 3. 散射振幅 M²

### 测试点
```
s = 500000 MeV² = 12.841 fm⁻²
t = -100000 MeV² = -2.568 fm⁻²
```

### Julia 结果
```
M²(uu_to_uu)       = 9.931266e+01 fm⁴
M²(us_to_us)       = 7.782093e+01 fm⁴
M²(ss_to_ss)       = 8.210192e+01 fm⁴
M²(ssbar_to_uubar) = 2.045641e+01 fm⁴
```

### Fortran 结果

**需要提取**：
1. 从散射振幅计算函数中提取 M² 的值
2. 使用相同的 s, t 值

**提取方法**：
- 找到计算 M² 的 Fortran 函数
- 添加 print 语句输出 M²(s, t)
- 对比相同 s, t 下的值

---

## 4. 截面 σ(s)

### 测试点

#### s = 200000 MeV² = 5.136 fm⁻²
```
σ(uu_to_uu)       = 2.431901e-02 fm²
σ(us_to_us)       = 2.148974e-02 fm²
σ(ss_to_ss)       = 2.278570e-02 fm²
σ(ssbar_to_uubar) = 5.366191e-02 fm²
```

#### s = 500000 MeV² = 12.841 fm⁻²
```
σ(uu_to_uu)       = 4.702861e-02 fm²
σ(us_to_us)       = 4.644946e-02 fm²
σ(ss_to_ss)       = 3.952626e-02 fm²
σ(ssbar_to_uubar) = 5.723492e-02 fm²
```

#### s = 1000000 MeV² = 25.682 fm⁻²
```
σ(uu_to_uu)       = 6.997242e-02 fm²
σ(us_to_us)       = 6.994222e-02 fm²
σ(ss_to_ss)       = 6.285315e-02 fm²
σ(ssbar_to_uubar) = 7.670286e-02 fm²
```

### Fortran 结果

**需要提取**：
1. 从 `tint` 函数（计算 σ(s)）中提取值
2. 使用相同的 s 值

**提取方法**：
- 在 `z1 tint.f90` 中添加 print 语句
- 输出 σ(s) 的值
- 对比相同 s 下的值

---

## 5. 平均散射率 ⟨Γ⟩

### Julia 结果
```
⟨Γ⟩(uu_to_uu)       = 5.614362e-02 fm⁻¹
⟨Γ⟩(ud_to_ud)       = 5.628311e-02 fm⁻¹
⟨Γ⟩(us_to_us)       = 5.713927e-02 fm⁻¹
⟨Γ⟩(ss_to_ss)       = 5.155648e-02 fm⁻¹
⟨Γ⟩(ssbar_to_uubar) = 1.895418e-01 fm⁻¹
```

### Fortran 结果

**需要提取**：
1. 从 `averaged_rate` 函数中提取 w_ij 的值
2. 对比单个散射过程的散射率

**提取方法**：
- 在 `z2 averaged_rate.f90` 中添加 print 语句
- 输出每个散射过程的 w_ij
- 对比相同过程的值

---

## 6. 弛豫时间 τ

### Julia 结果（已知）
```
τ_u = 1.726 fm
τ_s = 2.100 fm
τ_u/τ_s = 0.822
```

### Fortran 结果（已知）
```
τ_u = 0.581 fm
τ_s = 0.592 fm
τ_u/τ_s = 0.981
```

### 对比
```
τ_u 比值 = 2.97
τ_s 比值 = 3.55
```

---

## 🔍 调查策略

### 步骤 1: 验证数密度 ✅
- **结果**: 一致（差异 < 1%）
- **结论**: 数密度不是差异来源

### 步骤 2: 对比 K 值 ⏳
- **需要**: 从 Fortran 提取 A, G, K 值
- **方法**: 修改 `coefficients_tmu.f90`，添加输出
- **目标**: 验证有效耦合常数是否一致

### 步骤 3: 对比 M² ⏳
- **需要**: 从 Fortran 提取 M²(s, t) 值
- **方法**: 修改散射振幅计算函数，添加输出
- **目标**: 验证散射振幅计算是否一致

### 步骤 4: 对比 σ(s) ⏳
- **需要**: 从 Fortran 提取 σ(s) 值
- **方法**: 修改 `z1 tint.f90`，添加输出
- **目标**: 验证截面计算是否一致

### 步骤 5: 对比 ⟨Γ⟩ ⏳
- **需要**: 从 Fortran 提取单个过程的散射率
- **方法**: 修改 `z2 averaged_rate.f90`，添加输出
- **目标**: 找出哪个过程有差异

---

## 📝 Fortran 代码修改建议

### 1. 修改 `coefficients_tmu.f90`

在计算 A, G, K 后添加：
```fortran
! 输出中间量
print*, "=== Intermediate Values ==="
print*, "A_u = ", A_u
print*, "A_s = ", A_s
print*, "G_u = ", G_u
print*, "G_s = ", G_s
print*, "K_pi = ", K_pi
print*, "K_K = ", K_K
print*, "K_sigma = ", K_sigma
print*, "K_kappa = ", K_kappa
```

### 2. 修改散射振幅计算函数

添加测试代码：
```fortran
! 测试点
s_test = 500000.0d0 / hc**2  ! 500 GeV²
t_test = -100000.0d0 / hc**2  ! -100 GeV²

! 计算并输出 M²
M2_uu = calculate_M2_uu_to_uu(s_test, t_test)
print*, "M2(uu_to_uu, s=500GeV2, t=-100GeV2) = ", M2_uu
```

### 3. 修改 `z1 tint.f90`

在计算 σ(s) 后添加：
```fortran
! 输出测试点的截面
if (abs(s - 5.136d0) < 0.01d0) then  ! s = 200 GeV²
    print*, "sigma(s=200GeV2) = ", sigma
end if
```

### 4. 修改 `z2 averaged_rate.f90`

在计算 w_ij 后添加：
```fortran
! 输出单个过程的散射率
print*, "w_ij for process ", qk_code, " = ", w_ij
```

---

## 🎯 预期发现

根据目前的分析，差异可能来自：

1. **K 值的计算**：
   - 如果 K 值有差异，会导致 M² 和 σ(s) 都有差异
   - 这会直接影响弛豫时间

2. **M² 的计算**：
   - 如果 M² 有系统性差异（如因子差异）
   - 会导致 σ(s) 和弛豫时间都有差异

3. **σ(s) 的积分**：
   - 如果 t 积分有差异
   - 会导致 σ(s) 和弛豫时间有差异

4. **⟨Γ⟩ 的计算**：
   - 如果动量积分或角度积分有差异
   - 会直接影响弛豫时间

**最可能的情况**：
- 差异来自 **K 值** 或 **M²** 的计算
- 因为这会导致所有后续量都有系统性差异
- 3 倍的差异暗示可能有某个因子（如 2π, 4π 等）的差异

---

## 📚 相关文件

### Julia
- `deep_dive_fortran_julia_comparison.jl`: 计算 Julia 的中间量
- `AverageScatteringRate.jl`: 平均散射率计算
- `TotalCrossSection.jl`: 截面计算
- `ScatteringAmplitude.jl`: 散射振幅计算
- `EffectiveCouplings.jl`: 有效耦合常数计算

### Fortran
- `coefficients_tmu.f90`: 计算 A, G, K
- `z1 tint.f90`: 计算 σ(s)
- `z2 averaged_rate.f90`: 计算平均散射率
- `z1 relax_time.f90`: 计算弛豫时间

---

*文档时间: 2026-01-26*
*状态: 已完成 Julia 中间量计算，等待 Fortran 中间量提取*
