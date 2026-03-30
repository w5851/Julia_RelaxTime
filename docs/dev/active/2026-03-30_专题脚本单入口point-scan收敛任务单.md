# 专题脚本单入口 point-scan 收敛任务单

## 1. 目标

- 逐步实现“一个专题一个扫描脚本”，并在同一脚本支持 `mode=point|scan`。
- 在不破坏现有复现路径前提下压缩脚本数量与用户学习成本。

## 2. 前置依赖

- [ ] 求解器多初解 + 约束选优已合入主线。
- [ ] 回归门禁已归位到 `tests/regression/`。

## 3. 试点策略

- [ ] 试点专题：`scripts/pnjl/run_tmu_scan.jl`。
- [ ] 在试点内提供：
  - `mode=point`：单点求解/导出。
  - `mode=scan`：连续扫描。
- [ ] 兼容期保留旧脚本壳层，并输出迁移提示。

## 4. 实施步骤

- [ ] P1 定义统一参数契约（point 与 scan 的共享字段 + 专属字段）。
- [ ] P2 改造脚本 CLI 参数解析与默认值策略。
- [ ] P3 更新 `docs/guides/scripts/README.md` 与 `run_script_catalog.md`。
- [ ] P4 回填 `classification_manifest.toml` 与处置状态。
- [ ] P5 观察兼容周期后再处理弃用脚本删除决策。

## 5. 验证

- [ ] 单点 smoke 与扫描 smoke 均通过。
- [ ] 回归测试覆盖 point/scan 关键路径。
- [ ] 前端/CI 对接同一参数契约，不做双份映射。

## 6. 对接关系

- 上游：
  - `docs/dev/active/2026-03-30_求解器优化与脚本固化主计划.md`
  - `docs/dev/active/2026-03-30_回归门禁归位与脚本角色收敛任务单.md`
- 对接：
  - `docs/dev/active/2026-03-30_run脚本分类清单与清理候选任务单.md`
  - `docs/dev/active/2026-03-30_稳定脚本前端化与CI手动触发任务单.md`

## 7. 状态

- [ ] 待启动（建议作为 PR-3 核心内容）
