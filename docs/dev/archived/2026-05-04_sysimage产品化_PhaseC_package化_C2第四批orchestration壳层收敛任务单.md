---
title: sysimage 产品化 Phase C：package 化 C2 第四批 orchestration 壳层收敛任务单
archived: true
original: docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C2第四批orchestration壳层收敛任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase C：package 化 C2 第四批 orchestration 壳层收敛任务单

更新日期：2026-05-04

当前状态：第四批已完成；`run_gap_transport_scan.jl` 中重复的单点执行编排已收敛为独立 orchestration helper。

> 目的：在前几批已经拆出 CLI / IO / provenance / phase-equilibrium helper 的基础上，继续压薄 `run_gap_transport_scan.jl` 主体，把重复的“单点执行 + bulk fallback + 行输出”流程收敛成单一 orchestration shell。

---

## 1. 目标

- [x] C2-10 抽离 `run_scan` 单点 orchestration helper
- [x] C2-11 消除 xi 连续模式与常规扫描模式中的重复主流程
- [x] C2-12 保持原输出契约、失败容错与 continuity 语义不变

---

## 2. 本批范围

### 2.1 包含

- [x] 新增 `scripts/relaxtime/gap_transport_scan_orchestration.jl`
- [x] 更新 `scripts/relaxtime/run_gap_transport_scan.jl`
- [x] 更新 `tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl`

### 2.2 不包含

- [x] 不改动底层物理求解实现
- [x] 不改动 CSV schema / sidecar schema
- [x] 不直接把 orchestration helper 提升到 `src/` 稳定 API

---

## 3. 实现结果

### 3.1 orchestration helper

- [x] 新增 `GapTransportScanOrchestration` 模块
- [x] 抽离：
  - `build_scan_runtime`
  - `execute_gap_transport_scan_point!`

收益：

- `run_scan` 中原先两套重复的单点执行链现在合并为同一 helper
- 主脚本更清晰地区分“扫描循环控制”和“单点执行壳层”
- 后续若继续收敛，可直接围绕 orchestration helper 评估是否有第二调用方

### 3.2 稳定 API 对齐

- [x] orchestration helper 内部已改为优先走 `Models.solve_transport_from_equilibrium`

意义：

- 当平衡态已由 phase/equilibrium helper 求出后，不再回到 `solve_gap_and_transport(...; equilibrium=...)` 这一更厚的路径
- `run_gap_transport_scan.jl` 的 transport 后处理更贴近当前 `Models` 暴露的稳定 API 边界

---

## 4. 验证

- [x] `julia --project=. -e 'include("tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl")'`
- [x] `julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_run_scan_subcommand_dispatch_smoke.jl"; include("tests/integration/runtests.jl")'`
- [x] `julia --project=. scripts/dev/check_script_entrypoints.jl`
- [x] `julia --project=. scripts/dev/check_models_entry_contract.jl`

---

## 5. 当前脚本边界状态

截至本批，`scripts/relaxtime/run_gap_transport_scan.jl` 的脚本层已拆出：

- [x] CLI 参数解析
- [x] IO / header / sidecar row 写入
- [x] provenance 汇总壳层
- [x] phase tracker / equilibrium dispatch helper
- [x] 单点 orchestration helper

仍保留在主脚本中的主要职责：

- [ ] 扫描循环控制
- [ ] resume / overwrite / progress / gc 节奏管理
- [ ] 顶层入口装配

---

## 6. 下一步建议

- [ ] C2 收尾：评估扫描循环控制是否值得再抽一层 scan plan / iterator helper
- [ ] C3：让 precompile capability / workload 对齐已稳定的 transport point API，而不是仅覆盖脚本入口

---

## 7. DoD

- [x] 重复的单点执行链已收敛为单一 helper
- [x] continuity / failure handling / 输出契约保持不变
- [x] transport 后处理路径已更贴近稳定 API 边界
