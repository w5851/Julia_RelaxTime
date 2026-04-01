---
title: PNJL求解器解耦 Wave-B 任务单
archived: true
original: docs/dev/active/2026-04-01_PNJL求解器解耦Wave-B任务单.md
archived_date: 2026-04-01
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# PNJL求解器解耦 Wave-B 任务单

> 执行主线说明：本任务单用于 Wave-B 的唯一执行主线（勾选、证据、验收）。
> 设计与实现参考：
> - `docs/superpowers/specs/2026-04-01-pnjl-solver-decoupling-wave-b-design.md`
> - `docs/superpowers/plans/2026-04-01-pnjl-solver-decoupling-wave-b-implementation-plan.md`
> - `docs/architecture/models_solver_contract.md`

## 1. 目标

- [x] 将 Wave-A 已冻结契约推进到扫描/单点共用执行路径。
- [x] 准备兼容入口 deprecation 节奏（不在本波次强删）。
- [x] 固化 fallback/selection 诊断稳定性回归。

## 2. 范围

### 2.1 本期范围

- [x] scan 侧候选治理接线收敛。
- [x] 兼容入口统一路由与迁移状态补齐。
- [x] fallback 原因标签规范化。

### 2.2 非范围

- [x] 不移除 `solve_fixed*` 兼容入口。
- [x] 不替换求解后端。

## 3. 任务分解

### B1：scan 侧治理并轨

- [x] 为 scan 治理并轨补失败回归测试。
- [x] 让 `TmuScan/TrhoScan` 走共享 governance 接口。
- [x] 验证代表点无非预期漂移。

### B2：兼容路由 deprecation-ready

- [x] 为统一入口/兼容入口等价性补失败 smoke。
- [x] 完善迁移状态与查询接口。
- [x] 保持外部签名兼容。
- [x] 新建/维护 Wave-B 迁移映射台账：`docs/dev/active/2026-04-01_求解器兼容层迁移映射表_Wave-B.md`。

### B3：fallback 诊断稳定性

- [x] 为 fallback reason 标签补失败回归。
- [x] 统一 reason tag 与输出顺序。
- [x] 验证诊断输出可稳定复现。

### B4：文档治理同步

- [x] 更新迁移映射状态与删除门槛。
- [x] 回填验证证据到本任务单。
- [x] 通过 docs/governance 检查。

## 4. 验收标准

- [x] scan 与单点使用统一治理口径。
- [x] compatibility 路由行为可追溯且保持非破坏。
- [x] fallback 诊断回归稳定。
- [x] smoke + docs/governance 通过。

## 5. 验证命令

- [x] `julia --project=. -e 'include("tests/integration/models/test_waveb_compat_routing_smoke.jl")'`
- [x] `julia --project=. -e 'include("tests/regression/pnjl/test_waveb_scan_governance_stability.jl")'`
- [x] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [x] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [x] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [x] `julia --project=. scripts/dev/check_active_docs_governance.jl`

## 6. DoD

- [x] 任务项与验收项全部勾选。
- [x] 关键验证命令可复现通过。
- [x] 迁移状态文档同步更新。

## 7. 执行记录

- [x] 2026-04-01：基于 Wave-A 完成状态创建 Wave-B 任务单。
- [x] 2026-04-01：补充 Wave-B 迁移台账路径与新增测试命令，降低执行歧义。
- [x] 2026-04-01（B1）：新增 `tests/regression/pnjl/test_waveb_scan_governance_stability.jl` 并先失败后修复；`ScanCommon.attempt_with_candidates` 统一通过 governance 选择并输出 `governance.selection=*` 诊断。
- [x] 2026-04-01（B2）：新增 `tests/integration/models/test_waveb_compat_routing_smoke.jl` 并先失败后修复；新增 `Models.solver_migration_status` 查询接口，保持 legacy `solve_fixed*` 兼容入口可调用。
- [x] 2026-04-01（B3）：在同一 regression 中新增 fallback reason 标签归一化失败用例并修复；`ScanCommon.join_messages` 统一输出顺序 `governance.selection -> fallback.reason`。
- [x] 2026-04-01（验证）：定向测试、unit/integration/regression smoke、docs governance 均通过（integration smoke 期间可见预期 fail-on-fallback 演示日志，但测试总集 `Pass 353/353`）。
