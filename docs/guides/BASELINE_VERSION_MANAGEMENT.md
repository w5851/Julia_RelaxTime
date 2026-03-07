# Baseline Version Management Guide

## 概述

本文档规范了回归测试基线（baseline）的版本管理流程，确保在代码演进过程中能够安全地更新数值基线，同时保持测试的可追溯性和可复现性。

## 基线文件命名规范

### 当前规范

所有基线文件遵循以下命名格式：

```
baseline_{module}_{test_type}_v{version}.csv
```

**示例：**
- `baseline_pnjl_gap_fixedpoints_v1.csv`
- `baseline_njl_gap_fixedpoints_v1.csv`
- `baseline_transport_fixedpoints_v1.csv`

### 命名组件说明

1. **module**: 模块名称（pnjl, njl, njl2, rpnjl, phase, relaxtime 等）
2. **test_type**: 测试类型（gap_fixedpoints, constraint_fixedpoints, scan_fixedpoints 等）
3. **version**: 基线版本号（v1, v2, v3...）

## 基线版本演进流程

### 何时需要更新基线版本

在以下情况下需要创建新版本的基线：

1. **物理模型改进**：修正了物理模型的 bug 或改进了数值精度
2. **数值方法升级**：更换了求解器、积分方法或优化算法
3. **参数调整**：修改了模型参数（如耦合常数、截断能量等）
4. **覆盖范围扩展**：增加了新的固定点或测试场景
5. **容差策略调整**：基于物理理解重新评估了容差标准

### 不需要更新基线的情况

以下变更不应改变基线：

1. 代码重构（不改变数值结果）
2. 性能优化（在容差范围内保持数值一致）
3. 文档更新
4. 测试框架改进（不涉及数值计算）

## v1 → v2 迁移标准流程

### 第 1 步：验证当前基线

在创建新基线之前，确保当前基线测试全部通过：

```bash
# 运行完整回归测试
julia --project=. -e 'using Pkg; Pkg.test()'

# 或运行特定模块的回归测试
julia --project=. tests/regression/pnjl/test_pnjl_gap_fixedpoint_regression.jl
```

### 第 2 步：导出新基线

使用对应的导出脚本生成新版本基线：

```bash
# 示例：导出新的 PNJL gap 固定点基线
julia --project=. scripts/dev/export_pnjl_gap_fixedpoint_baseline.jl \
  --output tests/baselines/pnjl/baseline_pnjl_gap_fixedpoints_v2.csv
```

### 第 3 步：对比差异

使用 diff 或专用工具对比 v1 和 v2：

```bash
# 文本 diff（适合小文件）
diff tests/baselines/pnjl/baseline_pnjl_gap_fixedpoints_v1.csv \
     tests/baselines/pnjl/baseline_pnjl_gap_fixedpoints_v2.csv

# 或使用 Julia 脚本进行数值对比
julia --project=. scripts/dev/compare_baselines.jl \
  --old tests/baselines/pnjl/baseline_pnjl_gap_fixedpoints_v1.csv \
  --new tests/baselines/pnjl/baseline_pnjl_gap_fixedpoints_v2.csv \
  --rtol 1e-6 --atol 1e-10
```

### 第 4 步：更新测试文件

修改回归测试文件以使用新版本基线：

```julia
# 修改前
const GAP_BASELINE_PATH = joinpath(
    PROJECT_ROOT, "tests", "baselines", "pnjl",
    "baseline_pnjl_gap_fixedpoints_v1.csv"
)

# 修改后
const GAP_BASELINE_PATH = joinpath(
    PROJECT_ROOT, "tests", "baselines", "pnjl",
    "baseline_pnjl_gap_fixedpoints_v2.csv"
)
```

### 第 5 步：验证新基线

运行回归测试确保新基线有效：

```bash
julia --project=. tests/regression/pnjl/test_pnjl_gap_fixedpoint_regression.jl
```

### 第 6 步：记录变更

在 `docs/dev/active/` 下创建变更日志：

```markdown
# 文件：docs/dev/active/YYYY-MM-DD_baseline_v2_migration.md

## PNJL Gap Fixedpoints Baseline v1 → v2

**日期**：2026-03-06

**变更原因**：修正了 PNJL 求解器中的数值精度问题

**影响范围**：12 个固定点中的 3 个点数值发生变化

**最大相对差异**：5.2e-7（仍在物理可接受范围内）

**验证结果**：
- 所有回归测试通过
- 物理趋势保持一致
- 与文献数据对比无异常
```

### 第 7 步：归档旧基线

将旧版本基线移动到归档目录：

```bash
mkdir -p tests/baselines/archived/pnjl/
mv tests/baselines/pnjl/baseline_pnjl_gap_fixedpoints_v1.csv \
   tests/baselines/archived/pnjl/
```

**注意**：不要删除旧基线，保留用于历史对照和回溯。

### 第 8 步：更新文档

更新相关文档中的基线版本引用：

- `docs/dev/active/2026-03-06_项目全面审查与改进计划.md`
- `tests/regression/pnjl/README.md`（如果存在）
- CI 配置文件中的注释

### 第 9 步：提交变更

```bash
git add tests/baselines/pnjl/baseline_pnjl_gap_fixedpoints_v2.csv
git add tests/baselines/archived/pnjl/
git add tests/regression/pnjl/test_pnjl_gap_fixedpoint_regression.jl
git add docs/dev/active/YYYY-MM-DD_baseline_v2_migration.md
git commit -m "test(regression): migrate PNJL gap baseline from v1 to v2

- Reason: improved numerical precision in PNJL solver
- Max relative difference: 5.2e-7
- All regression tests pass with new baseline
- Archived v1 baseline for historical reference"
```

## 基线文件格式规范

### CSV 格式要求

1. **第一行必须是标题行**，包含所有字段名
2. **使用逗号分隔**，字段内不应包含逗号
3. **数值精度**：至少保留 16 位有效数字（Julia 默认 Float64 精度）
4. **注释行**：以 `#` 开头的行将被测试代码忽略

### 标准字段

每个基线 CSV 应包含以下标准字段（根据测试类型可能有所不同）：

**PNJL Gap Fixedpoints 示例：**
```csv
T,mu,phi_u,phi_d,phi_s,Phi,Phibar,omega,converged
```

**Transport Fixedpoints 示例：**
```csv
T,mu,xi,eta,sigma,zeta,tau_qq,tau_qg,converged
```

### 布尔值编码

- 使用 `1` 表示 `true`
- 使用 `0` 表示 `false`

### NaN 值处理

- 对于未定义或不适用的数值字段，使用 `NaN`（不带引号）
- 示例：约束求解中 `fixed_mu` 模式下 `rho_target` 应为 `NaN`

## 基线完整性检查清单

在提交新基线之前，确保：

- [ ] 所有固定点都已正确计算并收敛
- [ ] 数值精度达到要求（至少 16 位有效数字）
- [ ] CSV 格式符合规范（标题行、分隔符、编码）
- [ ] 测试代码能正确加载并解析基线文件
- [ ] 回归测试全部通过
- [ ] 已记录变更原因和验证结果
- [ ] 已归档旧版本基线

## 基线审查流程

### Pull Request 检查项

当 PR 涉及基线变更时，审查者应验证：

1. **变更合理性**：是否有充分理由更新基线？
2. **差异分析**：数值变化是否在物理可接受范围内？
3. **测试覆盖**：新基线是否通过了所有回归测试？
4. **文档完整**：是否记录了变更日志和验证结果？
5. **归档处理**：旧基线是否已正确归档？

### 自动化检查（未来计划）

考虑实现以下自动化检查：

- 基线格式验证脚本
- 数值差异自动报告
- 回归测试强制通过 gate
- 基线版本一致性检查

## 示例：完整的基线更新工作流

### 场景：修正 PNJL 模型中的数值 bug

```bash
# 1. 确保当前状态干净
git status
git checkout -b fix/pnjl-numerical-precision

# 2. 修复代码
# ... 编辑 src/models/pnjl_physics/PNJLCore.jl ...

# 3. 运行当前回归测试（预期失败）
julia --project=. tests/regression/pnjl/test_pnjl_gap_fixedpoint_regression.jl
# 预期输出：部分测试失败，数值偏差超出容差

# 4. 导出新基线
julia --project=. scripts/dev/export_pnjl_gap_fixedpoint_baseline.jl \
  --output tests/baselines/pnjl/baseline_pnjl_gap_fixedpoints_v2.csv

# 5. 对比差异
julia --project=. scripts/dev/compare_baselines.jl \
  --old tests/baselines/pnjl/baseline_pnjl_gap_fixedpoints_v1.csv \
  --new tests/baselines/pnjl/baseline_pnjl_gap_fixedpoints_v2.csv \
  --rtol 1e-6 --report baseline_diff_report.txt

# 6. 审查差异报告
cat baseline_diff_report.txt

# 7. 更新测试文件
sed -i 's/v1.csv/v2.csv/g' tests/regression/pnjl/test_pnjl_gap_fixedpoint_regression.jl

# 8. 验证新基线
julia --project=. tests/regression/pnjl/test_pnjl_gap_fixedpoint_regression.jl
# 预期输出：所有测试通过

# 9. 归档旧基线
mkdir -p tests/baselines/archived/pnjl/
mv tests/baselines/pnjl/baseline_pnjl_gap_fixedpoints_v1.csv \
   tests/baselines/archived/pnjl/baseline_pnjl_gap_fixedpoints_v1_deprecated_2026-03-06.csv

# 10. 记录变更
cat > docs/dev/active/2026-03-06_pnjl_gap_baseline_v2.md << EOF
# PNJL Gap Baseline v1 → v2 Migration

**Date**: 2026-03-06
**Reason**: Fixed numerical precision bug in PNJLCore.jl line 156
**Max Diff**: 5.2e-7 (relative)
**Status**: All regression tests pass
EOF

# 11. 提交变更
git add src/models/pnjl_physics/PNJLCore.jl
git add tests/baselines/pnjl/baseline_pnjl_gap_fixedpoints_v2.csv
git add tests/baselines/archived/pnjl/
git add tests/regression/pnjl/test_pnjl_gap_fixedpoint_regression.jl
git add docs/dev/active/2026-03-06_pnjl_gap_baseline_v2.md
git commit -m "fix(pnjl): improve numerical precision in gap solver

- Fixed floating-point accumulation error in PNJLCore.jl
- Updated baseline from v1 to v2 (max diff: 5.2e-7)
- All regression tests pass with new baseline
- Archived v1 baseline for historical reference"

# 12. 推送并创建 PR
git push origin fix/pnjl-numerical-precision
```

## 常见问题

### Q1: 如果我只是想试验性地测试新参数，是否需要创建 v2 基线？

**A**: 不需要。试验性测试应该：
1. 在本地生成临时基线（放在 `/tmp` 或 `tests/scratch/`）
2. 使用专门的测试脚本而不是修改回归测试
3. 只有在确定新参数成为默认配置后才更新正式基线

### Q2: 多个 PR 同时需要更新同一个基线怎么办？

**A**: 遵循以下原则：
1. **串行处理**：第一个 PR 更新到 v2，第二个 PR 更新到 v3
2. **协调合并**：如果可能，合并多个变更到一个 PR
3. **冲突解决**：使用最新的数值结果作为最终基线

### Q3: 如何处理向后兼容性测试？

**A**: 保留归档的旧基线：
```julia
# 可选：测试向后兼容性
@testset "Backward compatibility with v1 baseline" begin
    v1_baseline = load_baseline("archived/baseline_v1.csv")
    current_results = compute_results()
    # 验证趋势一致性而非精确数值
    @test correlation(v1_baseline, current_results) > 0.999
end
```

### Q4: 基线文件应该被 git 跟踪吗？

**A**: 是的。基线文件应该：
- ✅ 被 git 跟踪（提交到版本库）
- ✅ 参与代码审查
- ✅ 包含在 PR diff 中
- ❌ 不应该过大（如果单个基线 > 1MB，考虑拆分或压缩）

## 工具和脚本

### 当前可用的工具

1. **导出脚本**：`scripts/dev/export_*_baseline.jl`
2. **回归测试**：`tests/regression/**/test_*.jl`
3. **CI 集成**：`.github/workflows/nightly-full-regression.yml`

### 未来计划的工具

1. `scripts/dev/compare_baselines.jl` - 基线差异分析
2. `scripts/dev/validate_baseline_format.jl` - 格式验证
3. `scripts/dev/migrate_baseline.jl` - 自动化迁移助手

## 参考文档

- [项目全面审查与改进计划](../active/2026-03-06_项目全面审查与改进计划.md)
- [回归测试框架](../../../tests/regression/README.md)
- [CI 配置指南](.github/workflows/README.md)

---

**最后更新**：2026-03-06
**维护者**：Julia_RelaxTime 开发团队
