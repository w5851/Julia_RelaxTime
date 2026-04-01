# PNJL求解器解耦 Wave-C 任务单

> 执行主线说明：本任务单用于 Wave-C 的唯一执行主线（勾选、证据、验收）。
> 设计与实现参考：
> - `docs/superpowers/specs/2026-04-01-pnjl-solver-decoupling-wave-c-design.md`
> - `docs/superpowers/plans/2026-04-01-pnjl-solver-decoupling-wave-c-implementation-plan.md`
> - `docs/dev/archived/2026-03-31_PNJL求解器解耦框架约束与分层草案.md`

## 1. 目标

- [x] 将 Wave-B 统一治理能力延伸到扫描/工作流脚本主入口。
- [x] 收敛模型特化扫描 SOP 到 unified model-driven 路径。
- [x] 固化 Wave-C 迁移状态与稳定性回归基线。

## 2. 范围

### 2.1 本期范围

- [x] script/workflow 侧统一路由并轨。
- [x] 兼容适配层 migration 状态查询与台账同步。
- [x] model-driven 扫描输出稳定性回归。

### 2.2 非范围

- [ ] 不更换求解后端。
- [ ] 不在本波次执行兼容入口强删。

## 3. 任务分解

### C1：脚本与工作流统一路由

- [x] 为 script 路径与 unified 路径等价性补失败 smoke。
- [x] 将 scan 脚本入口收敛到统一 model-driven API。
- [x] 保持 CLI/外部调用签名兼容。

### C2：SOP 分叉收敛与迁移治理

- [x] 为兼容适配路由与迁移状态查询补失败测试。
- [x] 统一旧 SOP 为 compat adapter（unified-first）。
- [x] 更新/维护 Wave-C 扫描与工作流迁移映射台账：`docs/dev/active/2026-04-01_扫描与工作流兼容层迁移映射表_Wave-C.md`。

### C3：model-driven 稳定性回归

- [x] 为 representative points 补失败回归。
- [x] 统一诊断与输出顺序，保证断言稳定。
- [x] 验证回归可重复通过。

### C4：文档治理同步

- [x] 回填 Wave-C 任务单证据与勾选状态。
- [x] 同步迁移状态、删除门槛与风险备注。
- [x] 通过 docs/governance 检查。

## 4. 验收标准

- [x] 扫描/工作流新调用方默认走 unified model-driven 路径。
- [x] legacy SOP 保持 non-breaking 且行为可追溯。
- [x] model-driven 扫描回归稳定可复现。
- [x] smoke + docs/governance 通过。

## 5. 验证命令

- [x] `julia --project=. -e 'include("tests/integration/models/test_wavec_scan_workflow_parity_smoke.jl")'`
- [x] `julia --project=. -e 'include("tests/regression/pnjl/test_wavec_model_driven_scan_stability.jl")'`
- [x] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [x] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [x] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [x] `julia --project=. scripts/dev/check_active_docs_governance.jl`

## 5.1 默认零上下文 agent 执行约束

- [ ] 仅以下三份文档作为 Wave-C 执行主线：
  - `docs/dev/active/2026-04-01_PNJL求解器解耦Wave-C任务单.md`
  - `docs/superpowers/specs/2026-04-01-pnjl-solver-decoupling-wave-c-design.md`
  - `docs/superpowers/plans/2026-04-01-pnjl-solver-decoupling-wave-c-implementation-plan.md`
- [ ] 参考台账：`docs/dev/active/2026-04-01_扫描与工作流兼容层迁移映射表_Wave-C.md`
- [ ] 执行顺序必须严格按 C1 -> C2 -> C3（文档同步与验证在 C4/C5）。
- [ ] 每个子步采用 TDD（先失败测试，再最小实现，再回归）。
- [ ] 每个子步必须运行：本步定向测试 + unit/integration/regression smoke + docs governance checks。
- [ ] 不执行 Wave-C 之外的代码级改动；兼容路径只做 deprecation-ready，不做破坏性删除。

## 6. DoD

- [x] 任务项与验收项全部勾选。
- [x] 关键验证命令可复现通过。
- [x] Wave-C 迁移台账与状态文档同步更新。

## 7. 执行记录

- [x] 2026-04-01：创建 Wave-C 任务单草案，承接 archived 分层草案（Wave-C：扫描收敛阶段）与 Wave-B 已完成状态。
- [x] 2026-04-01：完成 C1→C2→C3（TDD）：新增 `tests/integration/models/test_wavec_scan_workflow_parity_smoke.jl` 与 `tests/regression/pnjl/test_wavec_model_driven_scan_stability.jl`，脚本入口 `run_tmu_scan/run_dense_trho_scan/run_adaptive_trho_scan` 收敛到 model-driven 参数化路由（含 `model_kind` 兼容参数）并保持 non-breaking；新增 `Models.scan_workflow_migration_map/status` 迁移查询面并同步台账。
- [x] 2026-04-01：完成 C4/C5 验证收口：定向测试 + unit/integration/regression smoke + docs governance 全部通过。
