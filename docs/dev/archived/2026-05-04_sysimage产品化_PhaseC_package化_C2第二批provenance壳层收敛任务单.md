---
title: sysimage 产品化 Phase C：package 化 C2 第二批 provenance 壳层收敛任务单
archived: true
original: docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C2第二批provenance壳层收敛任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase C：package 化 C2 第二批 provenance 壳层收敛任务单

更新日期：2026-05-04

当前状态：第二批已完成；`run_gap_transport_scan.jl` 尾部的 provenance/effective_config/summary/artifacts 汇总已抽离为独立 helper。

> 目的：延续 C2 首批的脚本壳层收敛，把 `run_gap_transport_scan.jl` 尾部的 sidecar 汇总逻辑从主执行体中移走，使主脚本进一步向 thin-wrapper + execution body 分离演进。

---

## 1. 目标

- [x] C2-4 抽离 provenance 汇总壳层
- [x] C2-5 为 helper 补最小 contract 测试
- [x] C2-6 保持原输出契约与入口不变

---

## 2. 本批范围

### 2.1 包含

- [x] 新增 `scripts/relaxtime/gap_transport_scan_provenance.jl`
- [x] 更新 `scripts/relaxtime/run_gap_transport_scan.jl`
- [x] 更新 `tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl`

### 2.2 不包含

- [x] 不改动 phase tracker / equilibrium / transport 主执行逻辑
- [x] 不改动 CSV schema、artifact 集合语义、run_manifest 语义
- [x] 不改动 `run_scan` 分发链

---

## 3. 实现结果

### 3.1 provenance helper

- [x] 新增 `GapTransportScanProvenance` 模块
- [x] 承载：
  - `build_effective_config`
  - `build_summary`
  - `collect_artifacts`
  - `write_scan_sidecars`

收益：

- `run_gap_transport_scan.jl` 不再直接展开大段 `effective_config / summary / artifacts` 构造
- sidecar 写入现在变成一个清晰的壳层调用点

### 3.2 单测补充

- [x] 在现有 unit 测试中补充 provenance shell contract
- [x] 校验：
  - `tau_interpolation_mode` / `integration_mode` 的字符串化口径
  - `summary` 计数汇总
  - `artifacts` 收集顺序

---

## 4. 验证

- [x] `julia --project=. -e 'include("tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl")'`
- [x] `julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_run_scan_subcommand_dispatch_smoke.jl"; include("tests/integration/runtests.jl")'`
- [x] `julia --project=. scripts/dev/check_script_entrypoints.jl`
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`

---

## 5. 当前脚本边界状态

截至本批，`scripts/relaxtime/run_gap_transport_scan.jl` 的脚本层已拆出：

- [x] CLI 参数解析
- [x] IO / header / sidecar row 写入
- [x] provenance 汇总壳层

仍保留在主脚本中的主要职责：

- [ ] phase tracker / hint
- [ ] equilibrium 调度
- [ ] transport scan 主执行体

---

## 6. 下一步建议

- [ ] C2-下一批：评估 `phase tracker + equilibrium dispatch` 是否值得形成稳定 helper / API
- [ ] C2-下一批：若继续收敛，优先抽“execution context / scan orchestration shell”，而不是直接拆细物理核
- [ ] C3：在稳定 API 边界更明确后，再调整 precompile workload 对齐

---

## 7. DoD

- [x] provenance 汇总壳层已抽离
- [x] helper 有最小 contract 测试
- [x] 原输出与入口契约保持不变
