---
title: sysimage 产品化 Phase C：package 化 C1 薄封装盘点任务单
archived: true
original: docs/dev/active/2026-05-04_sysimage产品化_PhaseC_package化_C1薄封装盘点任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase C：package 化 C1 薄封装盘点任务单

更新日期：2026-05-04

当前状态：已拉起并完成首轮盘点；稳定 CLI 顶层脚本厚度与 thin-wrapper 候选优先级已形成。

> 目的：在不触碰当前活跃物理主线的前提下，先识别哪些稳定 CLI 仍承载过厚的顶层逻辑，从而为后续 package 化时的入口收敛提供明确优先级。

---

## 1. 目标

- [x] C1 盘点当前稳定脚本中仍过厚的顶层逻辑
  - 验收：形成 thin-wrapper 候选列表

---

## 2. 盘点范围

稳定白名单入口：

- `scripts/pnjl/calculate_phase_structure.jl`
- `scripts/models/run_unified_scan.jl`
- `scripts/pnjl/run_conserved_charge_susceptibilities.jl`
- `scripts/relaxtime/run_relaxtime_orchestrator.jl`
- `scripts/relaxtime/run_gap_transport_scan.jl`
- `scripts/server/server_full.jl`

---

## 3. 盘点结果

### 3.1 脚本厚度概览

| 入口 | 行数 | 当前判断 | 说明 |
|---|---:|---|---|
| `scripts/relaxtime/run_gap_transport_scan.jl` | 1555 | 厚 | 兼有 CLI、扫描控制、输出、诊断、phase hint、求解调度等多类职责 |
| `scripts/pnjl/calculate_phase_structure.jl` | 518 | 中厚 | 虽然已调用 `Models.run_phase_pipeline`，但仍保留较多配置合并、manifest、CLI 规则 |
| `scripts/pnjl/run_conserved_charge_susceptibilities.jl` | 337 | 中厚 | 参数解析、scan 组织、CSV 输出、summary 打印仍在脚本层 |
| `scripts/models/run_unified_scan.jl` | 190 | 中等 | 已明显偏 thin-wrapper，但仍有一套手写参数解析与 phase workflow 转发 |
| `scripts/relaxtime/run_relaxtime_orchestrator.jl` | 97 | 薄 | 已基本是 package 化目标形态 |
| `scripts/server/server_full.jl` | 16 | 极薄 | 仅做入口转发，已接近最终形态 |

### 3.2 优先级建议

#### P0：最先收敛

- [x] `scripts/relaxtime/run_gap_transport_scan.jl`
  - 原因：体量最大，职责混杂最明显，也是 package 化收益最高的入口
  - 推荐拆分方向：
    - CLI 参数解析
    - 扫描配置 / provenance
    - equilibrium / phase hint / tracker
    - transport scan pipeline 主执行体

#### P1：第二批收敛

- [x] `scripts/pnjl/calculate_phase_structure.jl`
  - 原因：已经靠近 `Models`，但脚本层仍有较多 config/preset/manifest 逻辑
- [x] `scripts/pnjl/run_conserved_charge_susceptibilities.jl`
  - 原因：仍是“脚本内组织扫描 + 输出”的模式，尚未完全收敛为稳定 API 调用

#### P2：可择机处理

- [x] `scripts/models/run_unified_scan.jl`
  - 原因：当前已是相对薄的统一入口，可后续再视 CLI 契约抽象程度决定是否继续收敛

#### 当前不建议动

- [x] `scripts/relaxtime/run_relaxtime_orchestrator.jl`
- [x] `scripts/server/server_full.jl`
  - 原因：这两个入口已经基本满足 thin-wrapper 目标，不应为了 package 化而额外重构

---

## 4. 后续建议

### C2 主线建议

- [ ] 先从 `scripts/relaxtime/run_gap_transport_scan.jl` 启动最小收敛批次
- [ ] 目标不是一次性重写，而是先抽出稳定 API 壳层与脚本侧解析/输出边界
- [ ] 需要同步明确对应测试层，优先从 unit / integration 入口契约补证据

### C3 对 precompile 的启发

- [ ] 等稳定 API 边界更清晰后，再把 precompile workload 从“脚本偶然路径”进一步迁到“稳定 API 路径”

---

## 5. 验证依据

- [x] 读取稳定白名单脚本源码
- [x] 统计脚本行数
- [x] 盘点顶层 `function` / `include` / 入口结构

---

## 6. DoD

- [x] 稳定 CLI 厚度盘点完成
- [x] thin-wrapper 候选优先级完成
- [x] Phase C 下一步可直接从 `run_gap_transport_scan.jl` 启动
