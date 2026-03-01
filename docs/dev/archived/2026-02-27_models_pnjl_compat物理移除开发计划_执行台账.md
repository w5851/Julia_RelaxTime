---
title: models pnjl compat物理移除开发计划 执行台账
archived: true
original: docs/dev/active/2026-02-27_models_pnjl_compat物理移除开发计划_执行台账.md
archived_date: 2026-03-01
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# models pnjl compat物理移除开发计划 执行台账

更新日期：2026-02-27

---

## 0. 台账定位

- 本文档仅记录“执行事实”：变更点、命令、产物、结果。
- 本阶段开发任务定义统一维护在：
  - （未指定）

---

## 1. 记录规范（强约束）

- [x] 执行记录与开发任务分开保存：
  - 任务目标/范围/DoD 只写“开发任务单”；
  - 执行过程/命令/结果只写“执行台账”。
- [x] 追加记录时采用“直接追加”策略：
  - 不要求回读本台账历史上下文；
  - 仅按统一模板追加到文档末尾；
  - 每条记录必须自包含（含目标、命令、产物、结果）。
- [x] 每条记录必须可追溯到输出产物（`data/outputs/results/*`）。

---

## 2. 执行记录

- [ ] 批次号：Batch D4
  - 目标：完成 compat 能力迁移收口并物理删除目录，确保门禁全绿
  - 代码变更：迁移 compat core/solver/derivatives/analysis 与入口文件到 src/models/pnjl 主域；更新 ThermoFacade/EquilibriumFacade/PNJLMagneticModel 与 allowlist 路径；修复迁移后相对路径并删除 src/models/pnjl/compat 实体
  - 验证命令：
    - `UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`
    - `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD`
  - 产物：
    - `data/outputs/results/pnjl_prune_wave_snapshot_20260227_223601.txt`
    - `data/outputs/results/models_invokelatest_allowlist_delta_20260227_223601.txt`
    - `data/outputs/results/pnjl_scan_default_include_audit_20260227_223601.txt`
  - 结果：通过（smoke 845/845，三门禁 OK）
  - 主线映射（N1/N2/N3/N4）：D1/D2/D3/D4
