# 快速参考：Fortran vs Julia 弛豫时间对比

## 🎯 目标

找出 Fortran 和 Julia 弛豫时间 3 倍差异的根本原因。

---

## ✅ 已验证一致 (差异 < 2%)

| 物理量 | 状态 |
|--------|------|
| Polyakov 环 Φ | ✅ 0.26% |
| 有效质量 m | ✅ 1.75% |
| 单圈积分 A | ✅ 0.28% |
| 有效耦合 G | ✅ 2.0% |
| 数密度 ρ | ✅ 0.8% |

---

## 🔄 执行命令

### 1. 编译 Fortran
```powershell
cd d:\Desktop\fortran代码\输运系数\relaxtime_fortran
.\compile_gfortran.ps1
```

### 2. 运行 Fortran
```powershell
.\build\relaxtime.exe > results\debug_relaxtime.txt 2>&1
```

### 3. 运行 Julia 对比
```powershell
cd d:\Desktop\Julia_RelaxTime
julia scripts\relaxtime\extract_and_compare_relaxation_times.jl
```

---

## 📊 计算链

```
能隙方程 → (Φ, m) ✅
    ↓
A 函数 ✅
    ↓
G^f ✅
    ↓
K 系数 ❓
    ↓
M² ❓
    ↓
dσ/dt ❓
    ↓
σ(s) ❓
    ↓
w_ij ❓
    ↓
Γ ❓
    ↓
τ ❓
```

---

## 📁 关键文件

### Fortran
- `relaxtime_fortran/codes/relax time/z1 relax_time.f90` - 弛豫时间计算
- `relaxtime_fortran/codes/relax time/z2 averaged_rate.f90` - 平均散射率
- `relaxtime_fortran/results/debug_relaxtime.txt` - 调试输出

### Julia
- `Julia_RelaxTime/src/relaxtime/RelaxationTime.jl` - 弛豫时间模块
- `Julia_RelaxTime/src/relaxtime/AverageScatteringRate.jl` - 平均散射率
- `Julia_RelaxTime/scripts/relaxtime/extract_and_compare_relaxation_times.jl` - 对比脚本

### 文档
- `执行步骤.md` - 简化的执行步骤
- `NEXT_STEPS_RELAXTIME_COMPARISON.md` - 详细对比计划
- `COMPLETE_RESOLUTION.md` - 完整解决方案
- `FORTRAN_DEBUG_INSTRUCTIONS.md` - Fortran 调试指南

---

## 🎓 关键公式

### 弛豫时间
```
τ_u = 1 / Γ_u
Γ_u = Σ_j (n_j × w_uj)
```

### 平均散射率
```
w_ij = (N_c² / (2π⁴)) × ∫∫∫ d³p_i d³p_j dΩ × f_i × f_j × σ_ij × v_rel / (n_i × n_j)
```

### 总截面
```
σ(s) = ∫ dσ/dt dt
```

---

## 🔍 调试策略

### 如果 τ 一致 (< 10%)
✅ 问题解决! 之前的差异是参数或单位问题。

### 如果 τ 仍有差异 (> 50%)
逐步检查:
1. 对比所有 w_ij
2. 对比 σ(s) 在相同 s 值
3. 对比 M² 在相同动量
4. 对比 dσ/dt 在相同 s, t

---

*日期: 2026-01-26*
*状态: 准备运行对比*
