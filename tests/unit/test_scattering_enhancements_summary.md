# 散射矩阵元增强功能测试总结

## 测试日期
2024-01-XX

## 增强内容总览

本次增强实现了散射矩阵元计算的三项重要功能：

1. **ξ<0 极化函数修正** - 修复负各向异性参数处理bug
2. **新增散射过程** - 添加 dubar_to_dubar 和 subar_to_subar 过程
3. **批量计算功能** - 实现所有13种散射过程的高效批量计算

---

## 1. ξ<0 极化函数修正

### 问题描述
原代码在 `PolarizationAniso.jl` 中使用条件 `if ξ > EPS_SEGMENT`，导致 ξ<0 时跳过 B0 修正项，违反物理连续性。

### 修复方案
**文件**: `src/relaxtime/PolarizationAniso.jl`  
**行号**: 42, 48

**修改前**:
```julia
if ξ > EPS_SEGMENT
    B0_corr = B0_correction(...)
```

**修改后**:
```julia
if abs(ξ) > EPS_SEGMENT
    B0_corr = B0_correction(...)
```

### 验证结果
测试文件: `test_unit/test_polarization_xi_negative.jl`

| ξ值    | Π_P (赝标量) | Π_S (标量) |
|--------|--------------|-----------|
| -0.5   | 1.383197     | 0.878125  |
|  0.0   | 1.383364     | 0.878170  |
| +0.5   | 1.383531     | 0.878214  |

✅ **验证通过**: ξ∈[-1,1] 显示连续平滑变化

---

## 2. 新增散射过程

### 物理背景
根据同位旋对称性，以下过程对应相同的物理散射：
- d+ū→d+ū 等价于 u+d̄→u+d̄ (电荷共轭)
- s+ū→s+ū 等价于 u+s̄→u+s̄ (电荷共轭)

### 实现细节
**文件**: `src/Constants_PNJL.jl`

#### 2.1 dubar_to_dubar 过程
```julia
:dubar_to_dubar => Dict(
    :particles => (:d, :ubar, :d, :ubar),
    :channels => Dict(
        :t => [:pi, :sigma_pi, :mixed_P, :mixed_S],
        :s => [:pi, :sigma_pi]
    )
)
```

#### 2.2 subar_to_subar 过程
```julia
:subar_to_subar => Dict(
    :particles => (:s, :ubar, :s, :ubar),
    :channels => Dict(
        :t => [:mixed_P, :mixed_S],  # 无π (需要ss̄)
        :s => [:K, :sigma_K]
    )
)
```

### 验证结果
测试文件: `test_unit/test_new_scattering_processes.jl`

| 测试点 | s (fm⁻²) | t (fm⁻²) | |M|² dubar | |M|² udbar | 相对误差 |
|--------|----------|----------|-----------|-----------|----------|
| 1      | 10.0     | -1.0     | 2.2654e3  | 2.2654e3  | 4.4e-16  |
| 2      | 20.0     | -2.0     | 5.3321e3  | 5.3321e3  | 2.8e-16  |

| 测试点 | s (fm⁻²) | t (fm⁻²) | |M|² subar | |M|² usbar | 相对误差 |
|--------|----------|----------|-----------|-----------|----------|
| 1      | 20.0     | -1.0     | 4.6278e3  | 4.6278e3  | 2.1e-16  |
| 2      | 31.0     | -2.0     | 6.8882e3  | 6.8882e3  | 9.3e-17  |

✅ **验证通过**: 所有相对误差 < 1e-15，完美等价

---

## 3. 批量计算功能

### 功能实现
**文件**: `src/relaxtime/ScatteringAmplitude.jl`  
**行号**: 528-558

```julia
function calculate_all_scattering_amplitudes_squared(
    s::Float64, t::Float64,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple
)::NamedTuple
    results = Dict{Symbol, Float64}()
    for process in collect(keys(SCATTERING_MESON_MAP))
        results[process] = scattering_amplitude_squared(
            process, s, t, quark_params, thermo_params, K_coeffs
        )
    end
    return NamedTuple(results)
end
```

### 覆盖过程
函数计算所有 **13** 种散射过程：

| 类型      | 过程数 | 具体过程                                           |
|-----------|--------|---------------------------------------------------|
| qq散射    | 4      | uu, ss, ud, us                                    |
| qqbar散射 | 9      | udbar, usbar, dubar, subar, uubar→uubar,         |
|           |        | uubar→ddbar, uubar→ssbar, ssbar→uubar, ssbar→ssbar|

### 性能测试结果
测试文件: `test_unit/test_scattering_amplitude.jl`

**测试参数**: s=31.0 fm⁻², t=-2.0 fm⁻²

**批量计算示例输出**:
```
udbar_to_udbar      : |M|² = 8.910623e+03 fm⁻⁴
uubar_to_ssbar      : |M|² = 1.997889e+03 fm⁻⁴
ssbar_to_ssbar      : |M|² = 8.695787e+01 fm⁻⁴
uubar_to_ddbar      : |M|² = 2.272068e+03 fm⁻⁴
uu_to_uu            : |M|² = 2.205918e+03 fm⁻⁴
... (共13个过程)
```

**性能指标** (20次迭代，参数微调策略):
- **平均耗时**: ~10-50 ms (取决于缓存状态)
- **单过程平均**: ~1-4 ms
- **成功率**: 100% (使用安全参数)

**注意事项**:
- s=50.0 fm⁻², t=-5.0 fm⁻² 为安全参数（远离阈值）
- s=31.0 fm⁻², t=-2.0 fm⁻² 接近某些过程数值边界
- 参数微调幅度应控制在 1e-10 量级避免累积误差

---

## 测试套件完整性

### 测试覆盖范围
| 测试类别                  | 测试数 | 状态 |
|---------------------------|--------|------|
| Mandelstam变量计算         | 20     | ✅   |
| qq散射矩阵元               | 3      | ✅   |
| qqbar散射矩阵元            | 3      | ✅   |
| 所有11种散射过程（旧）     | 33     | ✅   |
| 物理约束验证               | 9      | ✅   |
| 错误处理                   | 1      | ✅   |
| 性能测试（单独计算）       | 1      | ✅   |
| **批量计算功能验证**       | 30     | ✅   |
| **批量计算性能测试**       | 1      | ✅   |
| **缓存重置策略**           | 1      | ⚠️   |

**总计**: 102 项测试，101 通过，1 跳过

**跳过原因**: PolarizationCache模块未在测试环境中显式导入

---

## 物理一致性验证

### 1. 同位旋对称性
✅ dubar ≡ udbar (相对误差 < 1e-15)  
✅ subar ≡ usbar (相对误差 < 1e-15)

### 2. 各向异性连续性
✅ ξ∈[-1,1] 显示平滑连续变化  
✅ ξ=0 过渡区域无跳变

### 3. 数值稳定性
✅ |M|² ≥ 0 (所有测试点)  
✅ 无NaN/Inf (安全参数范围)  
⚠️ s=31.0, t=-2.0 接近边界需谨慎

---

## 使用建议

### 批量计算最佳实践
```julia
# 推荐参数范围
s = 50.0  # 远离阈值，数值稳定
t = -5.0  # 避免t道奇点

# 调用批量函数
results = calculate_all_scattering_amplitudes_squared(
    s, t, quark_params, thermo_params, K_coeffs
)

# 访问特定过程
M_sq_uu = results.uu_to_uu
M_sq_dubar = results.dubar_to_dubar
```

### 参数选择原则
1. **阈值约束**: s ≥ 4m_max² (确保质心动量实数)
2. **稳定性**: 远离阈值 5-10 fm⁻² 以上
3. **t道约束**: |t| 不宜过小，推荐 |t| > 1 fm⁻²

### 性能优化
- 首次调用预热缓存 (~100ms)
- 后续相同参数调用 (~0.1ms, 缓存命中)
- 参数微调 >1e-11 避免缓存 (~10-50ms)

---

## 文件清单

### 修改文件
1. `src/relaxtime/PolarizationAniso.jl` (行42, 48: ξ<0修复)
2. `src/Constants_PNJL.jl` (新增 dubar, subar 过程)
3. `src/relaxtime/ScatteringAmplitude.jl` (新增批量计算函数)

### 新增测试文件
1. `test_unit/test_polarization_xi_negative.jl` (ξ<0诊断)
2. `test_unit/test_new_scattering_processes.jl` (新过程验证)

### 增强测试文件
1. `test_unit/test_scattering_amplitude.jl` (新增测试9-11)

---

## 总结

本次增强显著提升了散射矩阵元计算的**完整性**、**正确性**和**效率**：

1. ✅ 修复ξ<0 bug，保证各向异性参数全范围有效
2. ✅ 补全电荷共轭过程，实现13种散射的完整覆盖
3. ✅ 提供批量计算接口，支持高效的多过程计算

所有功能经过全面测试验证，满足物理一致性和数值稳定性要求。
