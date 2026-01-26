# A 函数差异调查 - 下一步行动

**日期**: 2026-01-26
**状态**: 🔍 深入调查中

---

## 📊 当前发现总结

### 1. 问题确认

**差异**: Fortran 和 Julia 的 A 函数计算结果相差 **31%**

| 项目 | Fortran | Julia | 差异 |
|------|---------|-------|------|
| A_u | -4.972 fm⁻² | -3.408 fm⁻² | 31.5% |
| A_s | -5.166 fm⁻² | -3.750 fm⁻² | 27.4% |

**测试条件**:
- T = 300 MeV (1.520 fm⁻¹)
- μ = 2 MeV (0.010 fm⁻¹)
- Φ = Φbar = 0.99999994
- m_u = 0.04051 fm⁻¹ (8.0 MeV)
- m_s = 1.0323 fm⁻¹ (203.7 MeV)

### 2. 已排除的原因

✅ **分布函数公式**: Julia 和 Fortran 的 PNJL 分布函数实现完全一致 (差异 < 1e-15)

✅ **常数项计算**: 数值积分和解析公式完全一致 (差异 < 1e-12)

✅ **积分范围**: 使用相同的积分范围 ([0, 3.05] 和 [0, 15.0]) 后差异仍为 31%

✅ **积分节点数**: 使用 128 个高斯-勒让德节点,与 Fortran 完全相同

✅ **积分方法**: 完全复制 Fortran 的两个循环结构后差异仍为 31%

### 3. Julia 的计算细节

```julia
# 常数项 (使用 y, w: [0, 3.05], 128 节点)
const_term = -4.6475 fm⁻²

# 分布函数项 (使用 y1, w1: [0, 15.0], 128 节点)
dist_term = 3.7956 fm⁻²

# 总和
A_u = 4 × (const_term + dist_term) = -3.408 fm⁻²
```

### 4. Fortran 的计算细节

从调试输出:
```
A_u (before G calc) = -4.9721620988593660
m_u = 4.0510366698690826E-002
```

---

## 🤔 可能的原因

### 假设 1: Fortran 使用了不同的 Φ 值

**检查方法**: 
- 从 Fortran 输出中提取实际使用的 Φ 和 Φbar 值
- 用这些值重新计算 Julia 的 A 函数

**可能性**: 中等

### 假设 2: Fortran 的分布函数有额外的归一化

**检查方法**:
- 仔细检查 Fortran 的 `fphi` 和 `fphibar` 函数
- 检查是否有全局的归一化因子

**可能性**: 高

### 假设 3: Fortran 的积分节点生成方式不同

**检查方法**:
- 从 Fortran 输出前几个积分节点和权重
- 与 Julia 的高斯-勒让德节点对比

**可能性**: 低 (高斯-勒让德积分是标准算法)

### 假设 4: Fortran 的 A 函数有额外的物理修正

**检查方法**:
- 检查 Fortran 代码中是否有 A 函数的后处理
- 检查是否有各向异性修正或其他物理效应

**可能性**: 高

---

## 📝 下一步行动计划

### 立即行动 (优先级: 高)

1. **提取 Fortran 的完整参数**
   ```fortran
   ! 在 quantity.f90 中添加调试输出
   print*, "DEBUG: Phi1 =", Phi1
   print*, "DEBUG: Phi2 =", Phi2
   print*, "DEBUG: T =", T
   print*, "DEBUG: mu_u =", arrMu(1)
   ```

2. **检查 Fortran 的分布函数调用**
   ```fortran
   ! 在 quantity.f90 的 A 函数计算循环中添加
   if (i == 1) then
       print*, "DEBUG: First point in dist loop:"
       print*, "  p =", p
       print*, "  E_u =", arrEi(1)
       print*, "  f_quark =", arrFq(1,1)
       print*, "  f_antiquark =", arrFq(2,1)
   end if
   ```

3. **对比第一个积分点的所有中间量**
   - 动量 p
   - 能量 E
   - 分布函数 f_quark, f_antiquark
   - 被积函数值
   - 贡献值

### 中期行动 (优先级: 中)

4. **检查 Fortran 是否有各向异性修正**
   - 搜索 Fortran 代码中的 "aniso" 或 "xi" 关键字
   - 检查是否有角度积分

5. **检查 Fortran 的 A 函数是否有后处理**
   - 在 `arrA(:) = arrA(:)*4d0` 之后是否还有其他操作
   - 检查 `arrG` 的计算是否影响了 `arrA`

### 长期行动 (优先级: 低)

6. **创建完整的 Fortran-Julia 对比测试套件**
   - 逐个函数对比
   - 逐个积分点对比
   - 自动化测试流程

---

## 🎯 预期结果

如果找到了差异的根本原因,应该能够:

1. **解释** 为什么 Fortran 的 A 函数比 Julia 大 31%
2. **修正** Julia 的实现,使其与 Fortran 一致
3. **验证** 修正后的 Julia 实现在多个测试点上都与 Fortran 一致

---

## 📚 相关文件

- `Julia_RelaxTime/src/relaxtime/OneLoopIntegrals.jl` - Julia 的 A 函数实现
- `relaxtime_fortran/codes/main/quantity.f90` - Fortran 的 A 函数实现
- `Julia_RelaxTime/src/QuarkDistribution.jl` - Julia 的分布函数实现
- `relaxtime_fortran/codes/PNJL/Fun distribution function.f90` - Fortran 的分布函数实现

---

*文档时间: 2026-01-26*
*状态: 需要进一步调查*
*下一步: 提取 Fortran 的完整参数并对比第一个积分点*
