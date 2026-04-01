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
