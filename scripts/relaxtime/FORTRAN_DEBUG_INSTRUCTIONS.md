# Fortran 调试输出指南

## 目标

在 Fortran 代码中添加调试输出,以便提取弛豫时间和中间计算结果。

---

## 第一步：添加弛豫时间输出

### 文件位置
`relaxtime_fortran/codes/relax time/z1 relax_time.f90`

### 修改内容

在 `subroutine relaxation_time` 的末尾,`arrTau` 赋值之前添加:

```fortran
! ============================================================================
! DEBUG: 输出弛豫时间
! ============================================================================
print*, "DEBUG: Relaxation times at T=", T*hc, "MeV, mu_B=", mu_B*hc, "MeV"
print*, "  tau_u  =", tau_l, "fm"
print*, "  tau_s  =", tau_s, "fm"
print*, "  tau_ub =", tau_lb, "fm"
print*, "  tau_sb =", tau_sb, "fm"
print*, "  Gamma_u  =", 1d0/tau_l, "fm^-1"
print*, "  Gamma_s  =", 1d0/tau_s, "fm^-1"
print*, "  Gamma_ub =", 1d0/tau_lb, "fm^-1"
print*, "  Gamma_sb =", 1d0/tau_sb, "fm^-1"
print*, ""
```

### 完整的修改后代码片段

```fortran
tau_l = 1d0 / (n_u*(w(2)+w(1))+n_ub*(w(6)+w(7)+w(9)+w(5))+n_s*w(3)+n_sb*w(8))
tau_s = 1d0 / (2d0*n_u*w(3)+2d0*n_ub*w(8)+n_s*w(4)+n_sb*(w(11)+2d0*w(10)))
tau_lb = 1d0 / (n_u*(w(6)+w(7)+w(9)+wa(5))+n_ub*(wa(1)+wa(2))+n_s*wa(6)+n_sb*wa(3))
tau_sb = 1d0 / (2d0*n_u*w(8)+2d0*n_ub*wa(3)+n_sb*wa(4)+n_s*(w(11)+2d0*w(10)))

! ============================================================================
! DEBUG: 输出弛豫时间
! ============================================================================
print*, "DEBUG: Relaxation times at T=", T*hc, "MeV, mu_B=", mu_B*hc, "MeV"
print*, "  tau_u  =", tau_l, "fm"
print*, "  tau_s  =", tau_s, "fm"
print*, "  tau_ub =", tau_lb, "fm"
print*, "  tau_sb =", tau_sb, "fm"
print*, "  Gamma_u  =", 1d0/tau_l, "fm^-1"
print*, "  Gamma_s  =", 1d0/tau_s, "fm^-1"
print*, "  Gamma_ub =", 1d0/tau_lb, "fm^-1"
print*, "  Gamma_sb =", 1d0/tau_sb, "fm^-1"
print*, ""

arrTau(1,:) = [tau_l, tau_l, tau_s]
arrTau(2,:) = [tau_lb, tau_lb, tau_sb]
end
```

---

## 第二步：添加平均散射率输出

### 文件位置
`relaxtime_fortran/codes/relax time/z2 averaged_rate.f90`

### 修改内容

在 `subroutine averaged_rate` 的末尾添加:

```fortran
! ============================================================================
! DEBUG: 输出平均散射率
! ============================================================================
! 注意：这会产生大量输出，建议只在需要时启用
! print*, "DEBUG: averaged_rate"
! print*, "  w_ij   =", w_ij, "fm^-1"
! print*, "  w_ij_n =", w_ij_n, "fm^-1"
! print*, ""
```

---

## 第三步：添加总截面输出

### 文件位置
`relaxtime_fortran/codes/relax time/z3 tint.f90` (或类似文件)

### 修改内容

在计算总截面的地方添加调试输出 (可选,用于详细调试):

```fortran
! DEBUG: 输出总截面
! print*, "DEBUG: sigma(s) at s=", s_cm, "fm^2"
! print*, "  sigma =", sig_ij, "fm^2"
```

---

## 第四步：重新编译

### Windows (使用 Intel Fortran)

```powershell
cd relaxtime_fortran
ifort /O2 /Qopenmp codes\*.f90 -o relaxtime.exe
```

### Linux/Mac (使用 gfortran)

```bash
cd relaxtime_fortran
gfortran -O2 -fopenmp codes/*.f90 -o relaxtime
```

或使用提供的编译脚本:

```bash
cd relaxtime_fortran
./compile_gfortran.sh
```

---

## 第五步：运行并保存输出

### Windows

```powershell
cd relaxtime_fortran
.\relaxtime.exe > results\debug_output.txt 2>&1
```

### Linux/Mac

```bash
cd relaxtime_fortran
./relaxtime > results/debug_output.txt 2>&1
```

---

## 第六步：提取调试信息

运行 Julia 脚本提取和对比:

```bash
cd Julia_RelaxTime
julia scripts/relaxtime/extract_and_compare_relaxation_times.jl
```

---

## 预期输出格式

Fortran 的调试输出应该类似:

```
DEBUG: Relaxation times at T=   29.595000000000002      MeV, mu_B=   78.920000000000002      MeV
  tau_u  =   1.2345678901234567E-002 fm
  tau_s  =   3.4567890123456789E-003 fm
  tau_ub =   1.2345678901234567E-002 fm
  tau_sb =   3.4567890123456789E-003 fm
  Gamma_u  =   81.000000000000000      fm^-1
  Gamma_s  =   289.35000000000002      fm^-1
  Gamma_ub =   81.000000000000000      fm^-1
  Gamma_sb =   289.35000000000002      fm^-1
```

---

## 注意事项

1. **单位**: Fortran 使用自然单位 (fm, fm⁻¹)
2. **输出位置**: 调试输出会打印到标准输出,需要重定向到文件
3. **性能**: 添加调试输出会略微降低性能,但对单点计算影响不大
4. **清理**: 调试完成后可以注释掉这些 print 语句

---

## 故障排除

### 问题 1: 编译错误

**症状**: `use` 语句找不到模块

**解决**: 确保所有依赖的模块文件都在编译路径中

### 问题 2: 运行时错误

**症状**: 程序崩溃或输出 NaN

**解决**: 检查输入参数是否正确,特别是温度和化学势

### 问题 3: 输出文件为空

**症状**: `debug_output.txt` 文件为空或很小

**解决**: 
- 检查程序是否正常运行完成
- 确保使用了输出重定向 `> results/debug_output.txt`
- 检查是否有错误输出 `2>&1`

---

*文档时间: 2026-01-26*
*状态: 准备添加 Fortran 调试输出*
