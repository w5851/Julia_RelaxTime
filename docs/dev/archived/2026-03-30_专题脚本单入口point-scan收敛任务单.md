---
title: 专题脚本单入口 point-scan 收敛任务单
archived: true
original: docs/dev/active/2026-03-30_专题脚本单入口point-scan收敛任务单.md
archived_date: 2026-03-31
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# 专题脚本单入口 point-scan 收敛任务单

## 1. 目标

- 逐步实现“一个专题一个扫描脚本”，并在同一脚本支持 `mode=point|scan`。
- 在不破坏现有复现路径前提下压缩脚本数量与用户学习成本。

约束决议（2026-03-30）：

- 目标范围为全部 `run_*.jl`，但执行方式是“按专题分波次压缩”。
- 本任务完成后，对应旧脚本立即删除，不保留兼容壳层。
- 回归门禁脚本不再保留在 `scripts/`，统一迁入 `tests/regression/`。

## 2. 前置依赖

- [x] 求解器多初解 + 约束选优已合入主线。
- [x] 回归门禁已归位到 `tests/regression/`。

## 3. 试点策略

- [x] 试点专题：`scripts/pnjl/run_tmu_scan.jl`。
- [x] 在试点内提供：
  - `mode=point`：单点求解/导出。
  - `mode=scan`：连续扫描。
- [x] 不保留兼容壳层；按波次提供迁移映射并执行旧脚本删除。

## 3.1 波次规划（按专题）

- [x] Wave-1（PNJL）：完成专题单入口压缩与旧脚本删除。
- [x] Wave-2（RelaxTime）：完成专题单入口压缩与旧脚本删除。
- [x] Wave-3（analysis/dev 等其余专题）：完成专题单入口压缩与旧脚本删除。

## 4. 实施步骤

- [x] P1 定义统一参数契约（point 与 scan 的共享字段 + 专属字段）。
- [x] P2 改造脚本 CLI 参数解析与默认值策略。
- [x] P3 更新 `docs/guides/scripts/README.md` 与 `run_script_catalog.md`。
- [x] P4 回填 `classification_manifest.toml` 与处置状态。
- [x] P5 对应波次脚本通过验证后，执行旧脚本即时删除。

## 5. 验证

- [x] 单点 smoke 与扫描 smoke 均通过。
- [x] 回归测试覆盖 point/scan 关键路径。
- [x] 前端/CI 对接同一参数契约，不做双份映射。
- [x] 对应专题的回归门禁断言已在 `tests/regression/` 落地，且 `scripts/` 无同类 regression 入口残留。

## 6. 对接关系

- 上游：
  - `docs/dev/active/2026-03-30_求解器优化与脚本固化主计划.md`
  - `docs/dev/active/2026-03-30_回归门禁归位与脚本角色收敛任务单.md`
- 对接：
  - `docs/dev/active/2026-03-30_run脚本分类清单与清理候选任务单.md`
  - `docs/dev/active/2026-03-30_稳定脚本前端化与CI手动触发任务单.md`

## 7. 状态

- [x] 已完成（PR-3 已落地：PNJL 试点单入口 `mode=point|scan` 与前后端契约接线完成）
