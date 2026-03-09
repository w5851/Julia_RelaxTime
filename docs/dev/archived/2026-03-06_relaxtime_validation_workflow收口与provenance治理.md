---
title: Relaxtime Validation Workflow 收口与 Provenance 治理
archived: true
original: docs/dev/active/2026-03-06_relaxtime_validation_workflow收口与provenance治理.md
archived_date: 2026-03-09
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Relaxtime Validation Workflow 收口与 Provenance 治理

## 背景

- 当前仓库中并存两条 `tau` 计算编排：
  - workflow 生产路径：`scripts/relaxtime/run_gap_transport_scan.jl` -> `TransportWorkflow.solve_gap_and_transport`
  - direct validation/scan 路径：`scripts/relaxtime/scan_relaxation_times_vs_T.jl` 与 `tests/validation/relaxtime/test_literature_digitized_tau_targets.jl`
- 两条路径在平衡态量（`Phi`, `m_s`, `n_s`）上基本一致，但在 `tau_s` 上出现系统分叉。
- 这导致“今天重新扫描得到的生产数据”与“validation/literature compare 结果”看起来互相矛盾。

## 已确认问题

- validation 之前没有复用生产 workflow，而是手工重建 `densities`、`sigma cache` 并直接调用 `RelaxationTime.relaxation_times`。
- `scan_relaxation_times_vs_T.jl` 的输出元数据不足，无法区分 `workflow` 与 `direct` 两条路径。
- `run_gap_transport_scan.jl` 在 fresh Julia 进程中缺少 `Models.jl` 入口加载，存在脚本上下文依赖。

## 本次立即修复

- 将 `tests/validation/relaxtime/test_literature_digitized_tau_targets.jl` 切换为基于 workflow 生产路径复算 literature targets。
- 修复 `scripts/relaxtime/run_gap_transport_scan.jl` 的 fresh-process 入口加载问题。
- 为 `run_gap_transport_scan.jl` 与 `scan_relaxation_times_vs_T.jl` 导出 CSV 增加 provenance 字段：
  - `git_commit`
  - `provenance.entrypoint`
  - `provenance.equilibrium_backend`
  - `provenance.tau_path`
  - `provenance.integration_mode`

## 后续建议

- 增加 workflow vs direct path 的一致性/差异性对照测试，明确哪些脚本允许分叉、哪些必须一致。
- 将 validation 分层：
  - production validation：面向 workflow 产线输出
  - kernel validation：面向底层 `RelaxationTime` 内核
- 为所有关键扫描 CSV 统一补充 provenance 最小集合，避免再次出现“同名文件、不同口径”的考古问题。

## 暂停前状态

- `tests/validation/relaxtime/test_literature_digitized_tau_targets.jl` 已切换到 workflow 生产路径，并在当前仓库状态下通过 `33/33` literature targets。
- `scripts/relaxtime/scan_relaxation_times_vs_T.jl` 已切换到 workflow 路径；重新生成的 `relaxation_times_vs_T_literature_compare.csv` 头部 provenance 已显示 `entrypoint=workflow`。
- `scripts/relaxtime/plot_tau_literature_comparison.py` 已能读取 scan provenance，并已生成新的对比图与 summary：
  - `data/outputs/figures/relaxtime/literature/tau_literature_comparison.png`
  - `data/outputs/figures/relaxtime/literature/tau_literature_comparison.pdf`
  - `data/outputs/results/relaxtime/scan/relaxation_times_vs_T_literature_compare_summary.csv`

## 下一步计划

- 补一组 workflow vs direct tau path 的显式对照测试，固定 `T=180/190/200 MeV, muB=0, xi=0` 等关键点，防止口径再次漂移却不自知。
- 评估是否保留 direct path：
  - 若仅作为内核诊断用途，则应从 production validation 与主绘图产线中剥离；
  - 若继续保留，则必须补充文档说明“与 workflow 不保证逐点一致”。
- 将 provenance 扩展到其余 relaxtime 输出脚本，至少覆盖 `git_commit`、`entrypoint`、`equilibrium_backend`、`tau_path`、`integration_mode`。
- 继续处理集成测试层面仍未完全收敛的问题，尤其是 clean process / 独立脚本入口的一致性验证。