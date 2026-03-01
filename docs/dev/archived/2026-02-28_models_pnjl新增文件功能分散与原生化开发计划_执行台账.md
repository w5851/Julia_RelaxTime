---
title: models pnjl新增文件功能分散与原生化开发计划 执行台账
archived: true
original: docs/dev/active/2026-02-28_models_pnjl新增文件功能分散与原生化开发计划_执行台账.md
archived_date: 2026-03-01
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# models pnjl新增文件功能分散与原生化开发计划 执行台账

更新日期：2026-02-28

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

- [ ] 批次号：Batch E1-E5
  - 目标：完成 pnjl 新增文件逐项审查并落地导数链路去桥接收口
  - 代码变更：冻结 src/models/pnjl 新增17文件清单并形成逐文件审查矩阵；新增 src/models/derivatives/ThermoDerivatives.jl；derivative_entrypoints 改为加载主域 ThermoDerivatives 模块；删除 ThermoDerivativesBridge.jl；同步 invokelatest allowlist
  - 验证命令：
    - `smoke='smoke'; julia --project=. tests/unit/runtests.jl`
    - `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD`
  - 产物：
    - `data/outputs/results/pnjl_prune_wave_snapshot_20260228_190255.txt`
    - `data/outputs/results/models_invokelatest_allowlist_delta_20260228_190255.txt`
    - `data/outputs/results/pnjl_scan_default_include_audit_20260228_190255.txt`
  - 结果：通过（smoke 845/845，四重门禁全绿）
  - 主线映射（N1/N2/N3/N4）：E1/E2/E3/E4/E5

- [ ] 批次号：Batch F1-F5
  - 目标：完成主域 scan 入口收口与 solver 空壳接口补齐，按顺序闭环 F1-F5
  - 代码变更：新增 src/models/scans/ScanEntrypoints.jl；新增 src/models/solver/{ConstraintModes,Solver}.jl 空壳接口；entrypoints/TransportWorkflow/MesonMassWorkflow 从 ScanEntrypointsBridge 切到 ScanEntrypoints；删除 ScanEntrypointsBridge.jl；同步 invokelatest allowlist
  - 验证命令：
    - `smoke='smoke'; julia --project=. tests/unit/runtests.jl`
    - `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD`
  - 产物：
    - `data/outputs/results/pnjl_prune_wave_snapshot_20260228_203222.txt`
    - `data/outputs/results/models_invokelatest_allowlist_delta_20260228_203222.txt`
    - `data/outputs/results/pnjl_scan_default_include_audit_20260228_203222.txt`
  - 结果：通过（smoke 845/845，四重门禁全绿）
  - 主线映射（N1/N2/N3/N4）：F1/F2/F3/F4/F5

- [ ] 批次号：Batch G1
  - 目标：将主域 solver 从模块引用壳推进为“真实类型 + 通用约束接口”
  - 代码变更：`src/models/solver/ConstraintModes.jl` 改为真实约束模式类型定义；`src/models/solver/Solver.jl` 改为 `solve_constraint(model, mode, T_fm; kwargs...)` 多重派发；`src/models/Models.jl` 接入 solver include 并导出 `solve_constraint`；移除 `Fixed*` 直接导出以消除 Main 名称冲突告警
  - 验证命令：
    - `$env:UNIT_PROFILE='smoke'; julia --project=. tests/unit/runtests.jl`
    - `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD`
  - 产物：
    - `data/outputs/results/pnjl_prune_wave_snapshot_20260228_224507.txt`
    - `data/outputs/results/models_invokelatest_allowlist_delta_20260228_224507.txt`
    - `data/outputs/results/pnjl_scan_default_include_audit_20260228_224507.txt`
  - 结果：通过（smoke 845/845，三门禁全绿）
  - 主线映射（N1/N2/N3/N4）：G1.1/G1.2/G1.3/G1.4/G1.5

- [ ] 批次号：Batch G2
  - 目标：将 TrhoScan 在 models backend 下的实际约束求解调用链收口到主域 `solve_constraint`
  - 代码变更：`src/models/scans/TrhoScan.jl` 新增 models 调用适配层（PNJL 约束类型 -> `Main.Models.Fixed*` 映射；`Main.Models.solve_constraint` 返回值 -> `SolverResult` 适配）；在 `_attempt_with_strategy` 与 `_solve_point` 中对 `solver_backend=:models` 切换到主域求解入口
  - 验证命令：
    - `julia --project=. tests/unit/pnjl/test_trho_scan_solver_backend_models_smoke.jl`
    - `$env:UNIT_PROFILE='smoke'; julia --project=. tests/unit/runtests.jl`
    - `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`
  - 产物：
    - `（本批未新增 prune-wave 产物，沿用 guard 命令输出）`
  - 结果：通过（trho models-backend 定向 smoke 通过；全局 smoke 845/845；migration/data-output guard 全绿）
  - 主线映射（N1/N2/N3/N4）：G2.1/G2.2/G2.3/G2.4

- [ ] 批次号：Batch G3
  - 目标：将 `TmuScan`（fixed-μ）与 `DualBranchScan` 在 models backend 下的实际约束求解调用链收口到主域 `solve_constraint`
  - 代码变更：`src/models/scans/TmuScan.jl` 新增 `FixedMu` 的 models 适配（`Main.Models.FixedMu` 映射 + `solve_constraint` 返回值适配为 `SolverResult`），并在 `_solve_point`/`_solve_point_with_seed_strategy` 切换到主域求解；`src/models/scans/DualBranchScan.jl` 新增 `solver_backend` 参数与同构 models 适配；两处均增加 fixed-μ models 关键词白名单透传（过滤 `iterations` 等不兼容参数）
  - 验证命令：
    - `julia --project=. tests/unit/pnjl/test_tmu_scan_solver_backend_models_smoke.jl`
    - `$env:UNIT_PROFILE='smoke'; julia --project=. tests/unit/runtests.jl`
    - `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`
  - 产物：
    - `（本批未新增 prune-wave 产物，沿用 guard 命令输出）`
  - 结果：通过（tmu models-backend 定向 smoke 通过；全局 smoke 845/845；migration/data-output guard 全绿）
  - 主线映射（N1/N2/N3/N4）：G3.1/G3.2/G3.3/G3.4/G3.5
