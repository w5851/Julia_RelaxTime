# PNJL 相图计算待办事项

> 创建日期：2025-12-30
> 状态：进行中

本文档记录 PNJL 相图计算相关的待完成任务，从 `PNJL_Solver_Refactoring_Notes.md` 和 `PNJL计算工作流.md` 中提取。

---

## 一、高优先级任务

### 1.1 Crossover 线计算

**背景**：在 T > T_CEP 区域，手征相变变为平滑的 crossover，需要定义和计算 crossover 线。

**定义方法**：

| 类型 | 方法 A（峰值法） | 方法 B（拐点法） |
|------|-----------------|-----------------|
| 手征 crossover | \|∂φ_u/∂T\| 峰值位置 | ∂²φ_u/∂T² = 0 的位置 |
| 退禁闭 crossover | \|∂Φ/∂T\| 峰值位置 | ∂²Φ/∂T² = 0 的位置 |

**技术基础**：
- 已有 `solve_with_derivatives()` 支持一阶和二阶导数计算
- `result.dx_dT[1]` = ∂φ_u/∂T，`result.dx_dT[4]` = ∂Φ/∂T
- `result.d2x_dT2[1]` = ∂²φ_u/∂T²，`result.d2x_dT2[4]` = ∂²Φ/∂T²

**算法流程**（固定 μ，扫描 T）：
1. 在 T ∈ [T_CEP, T_max] 范围内计算导数
2. 峰值法：找 |∂φ_u/∂T| 或 |∂Φ/∂T| 最大的 T
3. 拐点法：找 ∂²φ_u/∂T² 或 ∂²Φ/∂T² 符号变化的 T（二分法细化）

**实现计划**：
- [ ] 在 `PhaseTransition.jl` 中添加 `detect_crossover()` 函数
- [ ] 支持 `:chiral`（手征）和 `:deconfinement`（退禁闭）两种类型
- [ ] 支持 `:peak`（峰值法）和 `:inflection`（拐点法）两种方法
- [ ] 更新 `calculate_phase_structure.jl` 脚本
- [ ] 输出 `crossover.csv`

**输出格式**：
```
xi, mu_MeV, T_crossover_chiral_MeV, T_crossover_deconf_MeV
```

**注意事项**：
- 二阶导数计算成本较高（嵌套 AD）
- 二阶导数对数值噪声敏感，可能需要平滑处理
- 仅在 T > T_CEP 区域有效（crossover 区域）

### 1.2 T-μ 扫描脚本

**背景**：使用 `PhaseAwareContinuitySeed` 进行 T-μ 参数空间扫描。

**实现计划**：
- [ ] 创建 `scripts/pnjl/run_tmu_scan.jl`
- [ ] 支持命令行参数（ξ, T 范围, μ 范围）
- [ ] 自动加载 `boundary.csv` 用于相变感知
- [ ] 输出热力学量（P, s, ε, ρ, M_u, M_d, M_s）

### 1.3 热力学导数计算脚本

**背景**：计算质量导数和体粘滞系数。

**实现计划**：
- [ ] 创建 `scripts/pnjl/calculate_derivatives.jl`
- [ ] 计算 ∂M/∂T, ∂M/∂μ
- [ ] 计算体粘滞系数 ζ/s
- [ ] 输出 `derivatives.csv`, `bulk_viscosity.csv`

---

## 二、中优先级任务

### 2.1 GitHub Actions 自动化

**背景**：实现 PNJL 计算的自动化流水线。

**功能需求**：
- [ ] 手动触发（workflow_dispatch）
- [ ] 参数化输入（ξ, T 范围, μ/ρ 范围）
- [ ] 计算结果作为 Artifacts 上传
- [ ] 支持断点续算

**配置文件**：`.github/workflows/pnjl-pipeline.yml`

**工作流步骤**：
1. 安装 Julia 环境
2. 运行相结构计算
3. 运行参数空间扫描
4. 生成相图
5. 上传结果

### 2.2 完整测试套件

**背景**：新架构需要完整的测试覆盖。

**实现计划**：
- [ ] Core 模块单元测试
  - `Integrals.jl` 测试
  - `Thermodynamics.jl` 测试
- [ ] Solver 模块单元测试
  - `ConstraintModes.jl` 测试
  - `SeedStrategies.jl` 测试
  - `Conditions.jl` 测试
  - `ImplicitSolver.jl` 测试
- [ ] Derivatives 模块测试
- [ ] 集成测试
  - 扫描模块测试
  - 相变分析测试

---

## 三、低优先级任务

### 3.1 数值稳定性：Log-Sum-Exp 技巧

**背景**：在计算 Polyakov loop 修正的对数项时，指数函数可能导致数值溢出。

**参考实现**（Julia_test）：
```julia
function _log_polyakov_series(a, Phi1, Phi2)
    m = max.(0.0, a)
    m = max.(m, 2 .* a)
    m = max.(m, 3 .* a)
    term = exp.(-m) .+ 3.0 .* Phi1 .* exp.(a .- m) .+ 
           3.0 .* Phi2 .* exp.(2 .* a .- m) .+ exp.(3 .* a .- m)
    clamped = max.(term, POLYAKOV_EPS)
    return m .+ log.(clamped)
end
```

**当前状态**：使用 `safe_log` 处理，在极端参数下可能有问题。

**实现计划**：
- [ ] 评估当前实现在极端参数下的表现
- [ ] 如有必要，在 `Thermodynamics.jl` 中引入 log-sum-exp 技巧
- [ ] 添加数值稳定性测试

### 3.2 自适应动量截断

**背景**：当前使用固定动量上限（默认 10 fm⁻¹），在低温时可能浪费计算资源。

**Julia_test 实现**：
- 使用离散节点策略 [250, 500, 1000, 2000, 4000, 6000] MeV
- 通过 `ForwardDiff.value` 冻结参数，保证 AD 稳定性
- 低温快捷路径：T_mev ≤ (μ - m)/10 时直接返回 1.3 × pF

**权衡**：
- 固定上限：简单、AD 友好
- 自适应截断：更高效，但需要处理 AD 冻结带来的导数误差（~3-5%）

**实现计划**：
- [ ] 评估自适应截断的收益
- [ ] 如值得引入，在 `Integrals.jl` 中实现
- [ ] 确保与 ImplicitDifferentiation 兼容

### 3.3 诊断工具：Ω 组件分解

**背景**：便于调试和验证热力学势计算。

**参考实现**（Julia_test）：
```julia
function Omega_minimal_components(...)
    # 返回各组成部分
    return (
        chiral = Ω_chiral,
        polyakov = Ω_polyakov,
        vacuum = Ω_vacuum,
        thermal = Ω_thermal,
        total = Ω_total
    )
end
```

**实现计划**：
- [ ] 在 `Thermodynamics.jl` 中添加 `calculate_omega_components()` 函数
- [ ] 添加 `print_omega_components()` 辅助函数
- [ ] 用于调试和验证

---

## 四、可能的改进

### 4.1 代码重构

- [ ] 将 `calculate_phase_structure.jl` 中的重复代码迁移到模块
- [ ] 统一命令行参数解析

### 4.2 数据扩展

- [ ] 添加更多 ξ 值的数据（0.1, 0.3, 0.5）
- [ ] CEP 附近的细化扫描
- [ ] 高温区域的 crossover 数据

### 4.3 绘图改进

- [ ] 添加 3D 相图（T-μ-ξ）
- [ ] 添加热力学量等高线图
- [ ] 支持交互式绘图

---

## 五、依赖关系

```
Crossover 线计算
       ↓
T-μ 扫描脚本 ←── 依赖 boundary.csv + crossover.csv
       ↓
热力学导数计算
       ↓
GitHub Actions 自动化
```

---

## 六、相关文档

- `docs/dev/PNJL_Solver_Refactoring_Notes.md` - 重构笔记（已完成部分）
- `docs/notes/pnjl/PNJL计算工作流.md` - 工作流说明
- `docs/notes/pnjl/一阶相变方法对比.md` - 方法对比
- `docs/api/pnjl/PhaseTransition.md` - API 文档

---

*更新日期：2025-12-30*
