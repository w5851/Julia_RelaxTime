---
title: sysimage 产品化 Phase C：package 化 C3 首批 precompile/workload 对齐任务单
archived: true
original: docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C3首批precompile-workload对齐任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase C：package 化 C3 首批 precompile/workload 对齐任务单

更新日期：2026-05-04

当前状态：首批已完成；`scan/core/full` precompile profile 已新增基于稳定 transport API 的 capability。

> 目的：把 precompile coverage 从“更容易受脚本实现偶然路径影响”的层次，推进到更明确的稳定 API 边界，使 sysimage / warmup 更稳定命中当前 package 化后的 transport point 热路径。

---

## 1. 目标

- [x] C3-1 为 transport point 稳定 API 增加独立 precompile capability
- [x] C3-2 把 `scan/core/full` profile 对齐到该 capability
- [x] C3-3 维持现有 `scan_pipeline_cli` coverage，不回退已有扫描入口覆盖

---

## 2. 本批范围

### 2.1 包含

- [x] 更新 `src/models/precompile/registry.jl`
- [x] 更新 `scripts/dev/check_precompile_profile_coverage.jl`
- [x] 更新 `tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl`

### 2.2 不包含

- [x] 不移除现有 `scan_pipeline_cli` capability
- [x] 不直接对 sysimage trace budget 做结论性调整
- [x] 不把 script-private helper 直接提升为 `src/` 公共 API

---

## 3. 实现结果

### 3.1 新 capability

- [x] 新增 `:transport_point_api`

当前命中路径：

- [x] `Models.solve(...)`
- [x] `Models.solve_transport_from_equilibrium(...)`
- [x] `TransportWorkflow.TransportIntegrationConfig(...)`

设计判断：

- 该 capability 不再依赖 `run_gap_transport_scan.jl` 脚本入口
- 它覆盖的是当前更稳定的 transport point API 边界

### 3.2 profile 对齐

- [x] `scan` profile 已包含 `:transport_point_api`
- [x] `core / full / all` profile 已同步包含 `:transport_point_api`
- [x] `smoke / test` 保持原有轻量级策略不变

---

## 4. 验证

- [x] `julia --project=. scripts/dev/check_precompile_profile_coverage.jl`
- [x] `julia --project=. scripts/dev/check_models_entry_contract.jl`
- [x] `julia --project=. -e 'include("tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl")'`

---

## 5. 当前阶段判断

- [x] precompile/workload 已开始对齐稳定 transport API
- [x] 但 phase/equilibrium helper 仍保持 script-private，不应误判为已经完全 package 化

后续建议：

- [ ] 如需进一步压冷启动，可在后续 trace/benchmark 中比较 `transport_point_api` 加入前后的 trace residual
- [ ] 若出现第二个稳定调用方，再评估把更多 orchestration 语义提到 `src/models`

---

## 6. DoD

- [x] 稳定 transport API 已纳入 precompile capability
- [x] profile 覆盖矩阵已同步更新
- [x] 现有 scan pipeline coverage 未丢失
