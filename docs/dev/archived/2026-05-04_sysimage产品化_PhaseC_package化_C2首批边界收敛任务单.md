---
title: sysimage 产品化 Phase C：package 化 C2 首批边界收敛任务单
archived: true
original: docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C2首批边界收敛任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase C：package 化 C2 首批边界收敛任务单

更新日期：2026-05-04

当前状态：首批已完成；`run_gap_transport_scan.jl` 的 CLI 层与脚本 IO 层已从主执行体中抽离为独立 helper。

> 目的：以最小可合并批次启动 `run_gap_transport_scan.jl` 的 package 化收敛，先把“脚本壳层职责”和“扫描/求解执行体”分开，而不改动底层物理求解语义。

---

## 1. 目标

- [x] C2-1 抽离 CLI 参数解析层
- [x] C2-2 抽离脚本 IO / header / sidecar 层
- [x] C2-3 保持原入口与现有测试契约不变

---

## 2. 本批范围

### 2.1 包含

- [x] `scripts/relaxtime/gap_transport_scan_cli.jl`
- [x] `scripts/relaxtime/gap_transport_scan_io.jl`
- [x] `scripts/relaxtime/run_gap_transport_scan.jl` 改为组合 helper + 主执行体

### 2.2 不包含

- [x] 不改动 equilibrium / phase tracker / transport 核心求解逻辑
- [x] 不改动输出 schema
- [x] 不改动 `run_scan` 分发入口契约

---

## 3. 实现结果

### 3.1 CLI helper

- [x] 新增 `GapTransportScanCLI` 模块
- [x] 承载：
  - `ScanOptions`
  - `parse_args`
  - `print_usage`

收益：

- `run_gap_transport_scan.jl` 不再内嵌大段 CLI 参数解析
- 后续若继续 package 化，可直接把 CLI 和稳定 API 壳层进一步对接

### 3.2 IO helper

- [x] 新增 `GapTransportScanIO` 模块
- [x] 承载：
  - CSV header 写入
  - output header 兼容性检查
  - failed points sidecar 写入
  - 目录准备 / git commit 读取

收益：

- 主脚本中的“输出治理逻辑”与“扫描执行体”分离
- 后续若抽 `run_scan` 主执行壳层，边界更清晰

### 3.3 主脚本保留职责

当前 `scripts/relaxtime/run_gap_transport_scan.jl` 仍主要保留：

- [x] phase hint / tracker
- [x] equilibrium 调度
- [x] transport scan 主执行体
- [x] provenance sidecar 汇总调用

这符合“先拆脚本壳层，不碰求解语义”的本批边界。

---

## 4. 验证

- [x] `julia --project=. -e 'include("tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl")'`
- [x] `julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_run_scan_subcommand_dispatch_smoke.jl"; include("tests/integration/runtests.jl")'`
- [x] `julia --project=. scripts/dev/check_script_entrypoints.jl`
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`

---

## 5. 下一步建议

- [ ] C2-下一批：继续抽 `run_gap_transport_scan.jl` 的 provenance/effective_config 汇总壳层
- [ ] C2-下一批：评估 phase tracker / equilibrium 调度是否可形成更稳定的 package API 边界
- [ ] C3：待稳定 API 边界更清晰后，再调整 precompile workload 对齐新入口

---

## 6. DoD

- [x] CLI 层已从主脚本抽离
- [x] IO / sidecar 层已从主脚本抽离
- [x] 原入口与现有测试契约保持不变
