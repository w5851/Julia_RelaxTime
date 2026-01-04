# PNJL 相图计算待办事项

> 创建日期：2025-12-30
> 状态：进行中

本文档记录 PNJL 相图计算相关的待完成任务，从 `PNJL_Solver_Refactoring_Notes.md` 和 `PNJL计算工作流.md` 中提取。

---

## 一、高优先级任务

### 1.1 Crossover 线计算

**状态**：✅ 已实现

**背景**：在 T > T_CEP 区域，手征相变变为平滑的 crossover，需要定义和计算 crossover 线。

**定义方法**：

| 类型 | 方法 A（峰值法） | 方法 B（拐点法） |
|------|-----------------|-----------------|
| 手征 crossover | \|∂φ_u/∂T\| 峰值位置 | ∂²φ_u/∂T² = 0 的位置 |
| 退禁闭 crossover | \|∂Φ/∂T\| 峰值位置 | ∂²Φ/∂T² = 0 的位置 |

**实现位置**：`src/pnjl/analysis/PhaseTransition.jl`

**导出函数**：
- `detect_crossover(μ_fm, T_range; method, variable, xi, ...)` - 单点检测
- `scan_crossover_line(mu_range, T_range; ...)` - 扫描 crossover 线
- `CrossoverResult` - 结果结构体

**测试结果**（2025-12-30，修复后）：

| μ (MeV) | T_peak (MeV) | T_inflection (MeV) | 差异 (MeV) |
|---------|--------------|-------------------|------------|
| 0 | 200.15 | 200.16 | 0.003 |
| 50 | 198.27 | 198.27 | ~0 |
| 100 | 192.56 | 192.55 | 0.003 |
| 150 | 182.76 | 182.76 | ~0 |
| 200 | 168.67 | 168.67 | 0.005 |

**性能对比**：

| 测试 | 峰值法耗时 | 拐点法耗时 | 比例 |
|------|-----------|-----------|------|
| μ=0 MeV | ~34s | ~4s | 8.5x |
| μ=100 MeV | ~4s | ~4s | 1.0x |
| μ=200 MeV | ~2s | ~1s | 2.0x |

拐点法更快的原因：
- 峰值法使用黄金分割法，每次迭代需要 2 次导数计算
- 拐点法使用二分法，每次迭代只需要 1 次导数计算

**结论**：
1. 修复后两种方法在所有 μ 区域结果一致（差异 < 0.01 MeV）
2. **双峰结构**：手征凝聚 ∂φ_u/∂T 存在两个峰值：
   - 小峰（~147 MeV）：由退禁闭相变引起
   - 大峰（~200 MeV）：真正的手征 crossover
3. **多拐点结构**：∂²φ_u/∂T² 有多个符号变化：
   - 第一个 + → -（~150 MeV）：小峰的拐点
   - 最后一个 + → -（~200 MeV）：大峰的拐点（手征 crossover）
4. 算法改进（从高温向低温扫描）：
   - 峰值法：选择最大峰
   - 拐点法：选择第一个 - → +（高温方向的 + → -）
5. 性能：拐点法更快（~2-8x），推荐用于批量计算

**建议**：两种方法结果一致，推荐使用拐点法（更快）。

**测试脚本**：`scripts/pnjl/test_crossover_methods.jl`

**计算脚本集成**：`scripts/pnjl/calculate_phase_structure.jl`
- 添加 `--skip_crossover` 选项
- Step 5: 自动扫描 μ = 0 到 μ_CEP 的 crossover 线
- 输出文件：`crossover.csv`（格式：xi, mu_MeV, T_crossover_chiral_MeV, T_crossover_deconf_MeV）

**绘图脚本集成**：`scripts/pnjl/plot_phase_diagram.py`
- 添加 `--crossover` 和 `--no-crossover` 选项
- 手征 crossover：虚线（:）
- 退禁闭 crossover：点划线（-.）
- 与一阶相变线绘制在同一张图中

### 1.2 T-μ 扫描脚本

**状态**：✅ 已实现

**背景**：使用 `PhaseAwareContinuitySeed` 进行 T-μ 参数空间扫描。

**实现位置**：`scripts/pnjl/run_tmu_scan.jl`

**功能**：
- [x] 支持命令行参数（ξ, T 范围, μ 范围）
- [x] 自动加载 `boundary.csv` 用于相变感知
- [x] 输出热力学量（P, s, ε, ρ, M_u, M_d, M_s）
- [x] 断点续扫支持
- [x] 进度显示

**命令行参数**：
```
--xi=0.0          各向异性参数
--T_min=50        最低温度 (MeV)
--T_max=200       最高温度 (MeV)
--T_step=10       温度步长 (MeV)
--mu_min=0        最低化学势 (MeV)
--mu_max=400      最高化学势 (MeV)
--mu_step=10      化学势步长 (MeV)
--output=...      输出文件路径
--resume          断点续扫（默认启用）
--overwrite       覆盖已有文件
--no_phase_aware  禁用相变感知策略
--verbose         详细输出
```

**使用示例**：
```bash
# 基本扫描
julia scripts/pnjl/run_tmu_scan.jl

# 自定义参数范围
julia scripts/pnjl/run_tmu_scan.jl --xi=0.2 --T_min=50 --T_max=150 --mu_max=350

# 高分辨率扫描
julia scripts/pnjl/run_tmu_scan.jl --T_step=5 --mu_step=5 --verbose
```

**输出格式**：CSV 文件，包含列：
- T_MeV, mu_MeV, xi
- pressure_fm4, rho, entropy_fm3, energy_fm4
- phi_u, phi_d, phi_s, Phi1, Phi2
- M_u_MeV, M_d_MeV, M_s_MeV
- iterations, residual_norm, converged, message

### 1.3 热力学导数计算脚本

**状态**：✅ 已实现

**背景**：计算质量导数和体粘滞系数。

**实现位置**：`scripts/pnjl/calculate_derivatives.jl`

**功能**：
- [x] 创建 `scripts/pnjl/calculate_derivatives.jl`
- [x] 计算 ∂M/∂T, ∂M/∂μ
- [x] 计算体粘滞系数 ζ/s 相关量（v_n², ∂μ_B/∂T|_σ）
- [x] 输出 `derivatives.csv`, `bulk_viscosity.csv`

**命令行参数**：
```
--xi=0.0          各向异性参数
--T_min=100       最低温度 (MeV)
--T_max=200       最高温度 (MeV)
--T_step=10       温度步长 (MeV)
--mu_min=0        最低化学势 (MeV)
--mu_max=300      最高化学势 (MeV)
--mu_step=50      化学势步长 (MeV)
--output_dir=...  输出目录
--verbose         详细输出
```

**输出文件**：
1. `derivatives_xi{xi}.csv`：质量导数和热力学导数
   - T_MeV, mu_MeV, xi
   - M_u/d/s_MeV, dM_u/d/s_dT, dM_u/d/s_dmu
   - dP_dT, dP_dmu, dEpsilon_dT, dEpsilon_dmu, dn_dT, dn_dmu
   - pressure, energy, entropy, rho_norm

2. `bulk_viscosity_xi{xi}.csv`：体粘滞系数相关量
   - T_MeV, mu_MeV, xi
   - v_n_sq (等熵声速平方)
   - dmuB_dT_sigma (等比熵线斜率)
   - M_u/d/s_MeV, dM_u/d/s_dT, dM_u/d/s_dmuB
   - entropy, n_B

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

**状态**：✅ 已实现

**背景**：新架构需要完整的测试覆盖。

**实现位置**：`tests/unit/pnjl/`

**测试文件**：
- [x] `test_core_integrals.jl` - Core Integrals 模块测试（34 tests）
- [x] `test_core_thermodynamics.jl` - Core Thermodynamics 模块测试（37 tests）
- [x] `test_solver_constraint_modes.jl` - ConstraintModes 测试（27 tests）
- [x] `test_solver_seed_strategies.jl` - SeedStrategies 测试（66 tests）
- [x] `test_solver_conditions.jl` - Conditions 测试（28 tests）
- [x] `test_solver_implicit.jl` - ImplicitSolver 测试（64 tests）

**已有测试**（之前实现）：
- `test_phase_transition.jl` - 相变分析测试
- `test_scans.jl` - 扫描模块测试
- `test_thermo_derivatives.jl` - 热力学导数测试
- `test_bulk_viscosity.jl` - 体粘滞系数测试

**测试覆盖**：
- Core 模块：积分节点、真空项、热项、能量色散、热力学势、热力学量
- Solver 模块：约束模式、初值策略、条件函数、隐函数求解器
- 物理一致性：熵非负、质量正定、同位旋对称、P=-Ω

---

## 三、低优先级任务

### 3.1 数值稳定性：Log-Sum-Exp 技巧

**状态**：✅ 已完成

**背景**：在计算 Polyakov loop 修正的对数项时，指数函数可能导致数值溢出。

**实现**（2026-01-03）：

在 `src/pnjl/core/Integrals.jl` 中引入 Log-Sum-Exp 技巧：

```julia
@inline function log_polyakov_series_lse(a, Φ, Φ̄)
    m = a > 0 ? 3.0 * a : zero(a)
    term = exp(-m) + 3.0 * Φ * exp(a - m) + 3.0 * Φ̄ * exp(2.0 * a - m) + exp(3.0 * a - m)
    return m + log(max(term, POLYAKOV_EPS))
end
```

**性能测试结果**：

底层函数性能对比：
| 测试用例 | 旧版本(ns) | 新版本(ns) | 比值 |
|---------|-----------|-----------|------|
| 正常区域 | 142 | 217 | 1.52x |
| 相变区域 | 183 | 179 | 0.98x |
| **平均** | - | - | **1.45x** |

实际求解器性能（LSE 引入后）：
| 测试项目 | 耗时 | 收敛率 |
|---------|------|--------|
| 单点求解器（平均） | 9.45 ms | - |
| T-μ 扫描（49点） | 7.02 ms/点 | 100% |
| T-ρ 扫描（80点） | 54.32 ms/点 | 100% |

**数值稳定性测试**：

| 条件 | 结果 |
|------|------|
| 极低温 (T=10 MeV) | ✅ 收敛 |
| 大化学势 (μ=600 MeV) | ✅ 收敛 |
| 低温+大μ (T=30, μ=500) | ✅ 收敛 |
| 极高温 (T=500 MeV) | ✅ 收敛 |
| 极端条件 (T=20, μ=700) | ✅ 收敛 |

**结论**：
- 底层函数开销约 1.5x，但对实际求解器性能影响可忽略
- LSE 技巧成功解决了大化学势下的数值溢出问题
- 所有极端条件测试均通过，收敛率保持 100%

**实现计划**：
- [x] 评估当前实现在极端参数下的表现
- [x] 在 `Integrals.jl` 中引入 log-sum-exp 技巧
- [x] 添加数值稳定性测试

**相关文件**：
- 实现：`src/pnjl/core/Integrals.jl`
- 评估报告：`docs/dev/log_sum_exp_evaluation.md`
- 评估脚本：`scripts/pnjl/evaluate_log_sum_exp.jl`
- 性能测试：`scripts/pnjl/benchmark_lse_full.jl`

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

*更新日期：2026-01-03*
*1.2 实现日期：2025-12-31*
*1.3 实现日期：2025-12-31*
*2.2 实现日期：2026-01-02*
*3.1 实现日期：2026-01-03*
