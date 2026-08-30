---
title: Issue #130：phase-reference promotion gate v1 任务单
archived: true
original: docs/dev/active/2026-08-21_Issue130_phase-reference-promotion-gate.md
archived_date: 2026-08-30
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Issue #130：phase-reference promotion gate v1 任务单

状态：accepted；作者已接受 Issue #130 的 `strict_reference_v1`、
`derived_reference_v1` 和 `phase_surface_render_v1` 三层证据包。
promotion gate v1 为 `promotion_candidate`，但本任务单不执行 reference 写入、runtime
consumer 切换或 RS transport。

## 固定输入

- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`。
- numerical run：`32354095831`；同源 solver-free aggregate replay：`32451053476`。
- 三层 evidence：`docs/analysis/pnjl/phase_reference/issue130_phase_reference_layers_v1/`。
- promotion gate：`docs/analysis/pnjl/phase_reference/issue130_phase_reference_promotion_gate_v1/`。
- 审核 PR：`#248`，merge SHA：`09d1a6895100cef208c9108034d9f45b631158eb`。

## Gate 范围

`scripts/analysis/pnjl/audit_issue130_phase_reference_promotion.py` 以 solver-free
方式核验 calculation/source/replay provenance、三层完整性、strict unresolved 语义、
统一 161 个 xi、common-support interpolation、no-triangulation render、manifest/plot
hash 和作者审核记录。当前输出为 `promotion_candidate`；这表示可以准备独立的、
versioned phase-reference import PR，不表示已经完成晋升。

数值 CSV、PNG 和既有 `data/reference/pnjl` 保持不变。gate 明确记录
`reference_write=false`、`runtime_consumption=false` 和 `solver_called=false`（replay）。

## 收口与后续

1. promotion-gate 代码、审计包和本任务单已在 PR #250 合并，merge SHA 为
   `10231d83f8ed56b684f6acbd74a71ce85fdb47cf`。
2. 已转入独立的 versioned reference import task；该 PR 仅导入新的
   canonical sibling，并单独审计旧 reference 不变。
3. import PR 合并并完成 runtime consumer 评估后，才重新评估 RS transport。

## 非目标

- 不重跑 PNJL equilibrium solver，不触发新的 C0/C1/C2 或 Actions 数值运行。
- 不修改 Maxwell、endpoint policy、三态合同、几何容差或历史 evidence。
- 不因 `promotion_candidate` 自动写入 `data/reference` 或启动 transport。
