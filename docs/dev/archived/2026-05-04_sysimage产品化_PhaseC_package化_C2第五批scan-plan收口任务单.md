---
title: sysimage 产品化 Phase C：package 化 C2 第五批 scan-plan 收口任务单
archived: true
original: docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C2第五批scan-plan收口任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase C：package 化 C2 第五批 scan-plan 收口任务单

更新日期：2026-05-04

当前状态：第五批已完成；`run_gap_transport_scan.jl` 的扫描计划与循环控制已进一步收敛为独立 helper。

> 目的：在 orchestration helper 已抽离的基础上，继续压薄 `run_gap_transport_scan.jl` 顶层，明确“扫描计划生成 / 扫描循环控制 / 单点执行”的三层结构，让 C2 基本收口。

---

## 1. 目标

- [x] C2-13 抽离 scan plan 构造
- [x] C2-14 抽离扫描循环控制
- [x] C2-15 让主脚本更多退回到文件准备与顶层装配职责

---

## 2. 本批范围

### 2.1 包含

- [x] 新增 `scripts/relaxtime/gap_transport_scan_plan.jl`
- [x] 更新 `scripts/relaxtime/run_gap_transport_scan.jl`
- [x] 更新 `tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl`

### 2.2 不包含

- [x] 不改动物理求解、相态判据、CSV schema
- [x] 不把 scan plan helper 直接提升到 `src/` 稳定 API

---

## 3. 实现结果

### 3.1 scan plan helper

- [x] 新增 `GapTransportScanPlan` 模块
- [x] 抽离：
  - `build_scan_plan`
  - `execute_scan_plan!`

收敛内容：

- [x] `T_values / muB_values / xi_continuity_mode / total`
- [x] continuity 扫描与常规扫描的循环控制
- [x] progress / gc / skipped 计数节奏

### 3.2 主脚本当前职责

截至本批，`run_gap_transport_scan.jl` 顶层主要保留：

- [x] 文件与 sidecar 句柄准备
- [x] metadata/header 初始化
- [x] 调用 `build_scan_runtime / build_scan_plan / execute_scan_plan!`
- [x] 顶层 `main()` 入口

这意味着 C2 已基本达到“脚本更多承担入口与装配职责”的验收方向。

---

## 4. 验证

- [x] `julia --project=. -e 'include("tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl")'`
- [x] `julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_run_scan_subcommand_dispatch_smoke.jl"; include("tests/integration/runtests.jl")'`
- [x] `julia --project=. scripts/dev/check_script_entrypoints.jl`

---

## 5. 当前阶段判断

- [x] C2 已基本完成 `run_gap_transport_scan.jl` 的 script-private package 化收口
- [x] 当前剩余更像是“是否需要继续升格到 `src/` 公共 API”的决策问题，而不再是脚本内部结构混乱问题

后续建议：

- [ ] 若出现第二个稳定调用方，再评估把 scan plan / orchestration 的一部分迁移到 `src/models`
- [ ] Phase C 后续重点可转向 C3 的 precompile/workload 证据闭环

---

## 6. DoD

- [x] 扫描计划与循环控制已独立成 helper
- [x] 主脚本已基本退回装配职责
- [x] C2 包内收口目标基本达成
