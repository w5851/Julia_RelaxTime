---
title: models原生实现全仓替换legacy并移除执行台账
archived: true
original: docs/dev/active/2026-02-27_models原生实现全仓替换legacy并移除执行台账.md
archived_date: 2026-03-01
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# models原生实现全仓替换legacy并移除执行台账

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

- [ ] 批次号：Batch N1 + N2(部分)
  - 目标：完成入口层去 legacy 并启动调用方切换
  - 代码变更：新增 models 侧 scan/derivatives/core bridge；重构 entrypoints/derivative_entrypoints/workflows 路径；替换 smoke 相关 PNJL 单测中的 legacy_pnjl_module/Main.PNJL 调用；同步 models_invokelatest allowlist
  - 验证命令：
    - `UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`
    - `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`
  - 产物：
    - `（无）`
  - 结果：通过（N1 完成，N2 进行中）
  - 主线映射（N1/N2/N3/N4）：N1/N2

- [ ] 批次号：Batch N2
  - 目标：完成 tests/scripts 对 Models 原生 API 切换并清零 legacy 调用
  - 代码变更：批量替换 tests/** 与 scripts/** 中 Models.legacy_pnjl_module()/Main.PNJL 调用为 Models.pnjl_module()/PNJL 别名，保留非目标命名空间；同步通过 smoke 与 data-output 门禁
  - 验证命令：
    - `UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`
    - `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
  - 产物：
    - `（无）`
  - 结果：通过
  - 主线映射（N1/N2/N3/N4）：N2

- [ ] 批次号：Batch N3 + N4
  - 目标：完成 legacy 物理删除与全仓路径收敛
  - 代码变更：执行 git mv src/legacy_pnjl -> src/pnjl_compat 后，再下沉为 src/models/pnjl_compat；更新 src/models 与 tests 路径引用；移除 Models.legacy_pnjl_module 导出与入口别名；清零 src/tests/scripts 中 legacy_pnjl 与 Main.PNJL 直连
  - 验证命令：
    - `UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`
    - `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD`
    - `git ls-tree -r --name-only HEAD src/legacy_pnjl`
  - 产物：
    - `data/outputs/results/pnjl_prune_wave_snapshot_20260227_174155.txt`
    - `data/outputs/results/models_invokelatest_allowlist_delta_20260227_174155.txt`
    - `data/outputs/results/pnjl_scan_default_include_audit_20260227_174155.txt`
  - 结果：通过（src/legacy_pnjl 已物理删除）
  - 主线映射（N1/N2/N3/N4）：N3/N4

- [ ] 批次号：Batch N3 收尾（compat 下沉）
  - 目标：完成 DoD 第二项：主路径全部走 models，pnjl_compat 依赖下沉到 models
  - 代码变更：执行 git mv src/pnjl_compat -> src/models/pnjl_compat；收敛 src/models/** 与 tests/** 路径常量；修复 src/models/pnjl_compat/** 内部相对路径层级；同步 models_invokelatest_allowlist 新路径条目
  - 验证命令：
    - `UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`
    - `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD`
  - 产物：
    - `data/outputs/results/pnjl_prune_wave_snapshot_20260227_182912.txt`
    - `data/outputs/results/models_invokelatest_allowlist_delta_20260227_182912.txt`
    - `data/outputs/results/pnjl_scan_default_include_audit_20260227_182912.txt`
  - 结果：通过（主路径已收敛至 models；compat 目录已下沉）
  - 主线映射（N1/N2/N3/N4）：N3

- [ ] 批次号：Batch N3 架构对齐（pnjl 子域）
  - 目标：将 compat 目录并入 models/pnjl 子域并完成路径收敛
  - 代码变更：执行 git mv src/models/pnjl_compat -> src/models/pnjl/compat；更新 src/models/core|derivatives|scans|legacy|pnjl 与 tests 路径常量；修复 src/models/pnjl/compat/** 内部相对路径层级；更新 models_invokelatest_allowlist 到新路径
  - 验证命令：
    - `UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`
    - `julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`
    - `julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD`
  - 产物：
    - `data/outputs/results/pnjl_prune_wave_snapshot_20260227_192330.txt`
    - `data/outputs/results/models_invokelatest_allowlist_delta_20260227_192330.txt`
    - `data/outputs/results/pnjl_scan_default_include_audit_20260227_192330.txt`
  - 结果：通过（目录结构对齐 models/pnjl，回归与门禁全绿）
  - 主线映射（N1/N2/N3/N4）：N3
