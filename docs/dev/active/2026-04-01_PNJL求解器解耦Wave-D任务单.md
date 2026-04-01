# PNJL求解器解耦 Wave-D 任务单

> 执行主线说明：本任务单用于 Wave-D 的唯一执行主线（勾选、证据、验收）。
> 设计与实现参考：
> - `docs/superpowers/specs/2026-04-01-pnjl-solver-decoupling-wave-d-design.md`
> - `docs/superpowers/plans/2026-04-01-pnjl-solver-decoupling-wave-d-implementation-plan.md`
> - `docs/dev/archived/2026-03-31_PNJL求解器解耦框架约束与分层草案.md`

## 1. 目标

- [ ] 在阈值满足前提下完成 compat 入口清理（删除或硬弃用）。
- [ ] 保持 solver/scan/workflow 统一 model-driven 主路径唯一性。
- [ ] 固化 Wave-D 清理后治理可追溯性与稳定性回归基线。

## 2. 范围

### 2.1 本期范围

- [ ] solver + scan/workflow 兼容入口阈值化清理。
- [ ] migration 状态由 deprecation-ready 收敛到 removed/hard-deprecated。
- [ ] 清理后 deterministic 回归与 smoke 验证。
- [ ] 稳定公共入口变更同步更新 `docs/api/`。

### 2.2 非范围

- [ ] 不更换求解后端。
- [ ] 不新增无关模型功能范围。

## 3. 任务分解

### D1：失败测试先行（兼容清理契约）

- [ ] 为 legacy 入口清理后的期望行为补失败 integration 测试。
- [ ] 为清理后输出稳定性补失败 regression 测试。
- [ ] 验证失败原因与 Wave-D 目标一致。

### D2：阈值化清理执行（删除或硬弃用）

- [ ] `solve_fixed*` 及高风险 compat 路由默认先 hard-deprecate，阈值满足后再 remove。
- [ ] 保留 migration 状态查询可用并标注 removed/hard-deprecated。
- [ ] 确保统一主路径行为不变（unified-first）。

### D3：迁移台账与治理同步

- [ ] 更新 Wave-D 迁移映射台账：`docs/dev/active/2026-04-01_扫描与工作流兼容层迁移映射表_Wave-D.md`。
- [ ] 回填删除门槛满足证据与回退说明。
- [ ] 保持 docs 治理检查通过。

### D4：验证矩阵与收口

- [ ] 运行定向 Wave-D integration/regression。
- [ ] 运行 unit/integration/regression smoke。
- [ ] 运行 docs governance checks。
- [ ] 运行 `scripts/dev/check_unit_skip_policy.jl` 并记录结果。
- [ ] 回填任务单证据与勾选状态。

## 4. 验收标准

- [ ] compat 清理动作均有阈值证据与迁移映射可追溯。
- [ ] unified model-driven 主路径为唯一稳定入口。
- [ ] 清理后 regression/smoke 结果稳定可复现。
- [ ] docs/governance 通过。

## 5. 验证命令

- [ ] `julia --project=. -e 'include("tests/integration/models/test_waved_compat_cleanup_smoke.jl")'`
- [ ] `julia --project=. -e 'include("tests/regression/pnjl/test_waved_cleanup_stability.jl")'`
- [ ] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [ ] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [ ] `julia --project=. scripts/dev/check_active_docs_governance.jl`
- [ ] `julia --project=. scripts/dev/check_unit_skip_policy.jl`
- [ ] `docs/api/` 变更核对（如涉及稳定公共入口）

## 6. DoD

- [ ] 任务项与验收项全部勾选。
- [ ] 关键验证命令可复现通过。
- [ ] Wave-D 迁移台账与状态文档同步更新。

## 7. 执行记录

- [x] 2026-04-01：创建 Wave-D 三件套草案（active/spec/plan），基于 2026-03-31 分层草案与 Wave-C 收口现状，作为下一 PR 的执行主线。

## 8. 本轮执行证据（2026-04-01）

### D1：失败测试先行（兼容清理契约）

- [x] 为 legacy 入口清理后的期望行为补失败 integration 测试。
- [x] 为清理后输出稳定性补失败 regression 测试。
- [x] 验证失败原因与 Wave-D 目标一致。

证据（RED）：

- `julia --project=. -e 'include("tests/integration/models/test_waved_compat_cleanup_smoke.jl")'`
  - 失败点：`status.status` 仍为 `:active`；`solve_fixedmu_constraint` 尚未 hard-deprecated。
- `julia --project=. -e 'include("tests/regression/pnjl/test_waved_cleanup_stability.jl")'`
  - 失败点：`solve_fixedmu_constraint` 尚未抛 hard-deprecate 错误。

### D2：阈值化清理执行（删除或硬弃用）

- [x] `solve_fixed*` 及高风险 compat 路由默认先 hard-deprecate，阈值满足后再 remove。
- [x] 保留 migration 状态查询可用并标注 removed/hard-deprecated。
- [x] 确保统一主路径行为不变（unified-first）。

证据（GREEN）：

- `julia --project=. -e 'include("tests/integration/models/test_waved_compat_cleanup_smoke.jl")'` 通过（7/7）。
- `julia --project=. -e 'include("tests/regression/pnjl/test_waved_cleanup_stability.jl")'` 通过（6/6）。
- 回归兼容测试同步：
  - `julia --project=. -e 'include("tests/integration/models/test_waveb_compat_routing_smoke.jl")'` 通过（9/9）。
  - `julia --project=. -e 'include("tests/unit/models/test_solver_compat_markers.jl")'` 通过（16/16）。
  - `julia --project=. -e 'include("tests/unit/models/test_constraint_solver.jl")'` 通过（39/39）。

### D3：迁移台账与治理同步

- [x] 更新 Wave-D 迁移映射台账：`docs/dev/active/2026-04-01_扫描与工作流兼容层迁移映射表_Wave-D.md`。
- [x] 回填删除门槛满足证据与回退说明。
- [x] 保持 docs 治理检查通过。

证据：

- 迁移映射表四类 compat 路由状态已更新为 `hard_deprecated`。
- `docs/api/` 同步更新：
  - `docs/api/models/solver/ConstraintModes.md`
  - `docs/api/models/solver/Overview.md`
  - `docs/api/models/solver/README.md`
- 文档治理：
  - `julia --project=. scripts/dev/check_docs_consistency.jl` -> OK
  - `julia --project=. scripts/dev/check_active_docs_governance.jl` -> OK

### D4：验证矩阵与收口

- [x] 运行定向 Wave-D integration/regression。
- [x] 运行 unit/integration/regression smoke。
- [x] 运行 docs governance checks。
- [x] 运行 `scripts/dev/check_unit_skip_policy.jl` 并记录结果。
- [x] 回填任务单证据与勾选状态。

验证摘要：

- 定向：Wave-D integration/regression 均通过。
- smoke：
  - unit smoke 通过（781/781）。
  - integration smoke 通过（355/355）。
  - regression smoke 通过（462 pass，1 broken；broken 为既有可选夹具项）。
- 治理：
  - `check_unit_skip_policy` -> OK
  - `check_docs_consistency` -> OK
  - `check_active_docs_governance` -> OK
