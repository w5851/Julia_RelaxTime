---
title: Phase pipeline baseline v3：严格三交点 Maxwell 语义
archived: true
original: docs/dev/active/2026-08-16_phase_pipeline_baseline_v3.md
archived_date: 2026-08-16
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Phase pipeline baseline v3：严格三交点 Maxwell 语义

创建日期：2026-08-16

状态：本 PR 已完成本地 candidate 生成与 focused regression 验证，等待 GitHub Actions
验证后合并。旧 `baseline_phase_pipeline_v2.csv` 保留为历史 contract，不覆盖、不删除。

## 变更原因

公共 `PhaseCore.maxwell_construction` 当前默认使用
`candidate_policy=:unique_three_crossing_sign_change_v2`。该语义只接受恰有三个去重交点、
且具有非零宽度面积变号 bracket 的 Maxwell 候选；旧
`:unique_three_crossing_topology_v1` 仅用于历史 artifact 重放。

`baseline_phase_pipeline_v2.csv` 由 `be7545d9`（2026-07-19）引入，早于
`ceec2295`（PR #174，2026-08-05）和 `d1f2b52`（2026-08-13）的 Maxwell 语义收口。
因此本次创建 v3，而不是无审计地覆盖 v2。

## 影响范围

- baseline：新增 `tests/baselines/phase/baseline_phase_pipeline_v3.csv`。
- 测试：phase regression 改用 v3；断言容差保持 `rtol=1e-6`、`atol=1e-10`。
- 导出：导出脚本默认输出改为 v3。
- 代码：不修改求解器、Maxwell 算法、CEP 判定或 crossover 算法。

## v2 到 v3 差异

候选由 `origin/main@bf249323` 的导出脚本生成。两次独立导出字节级一致，SHA-256 为
`42470B014971C36D5BD419C8EAF2118F4090398D5AFEC73B004A39D9FF0B9E56`。

- `T=120 MeV` boundary 行保持不变。
- `T=125 MeV` boundary：
  - `mu_transition_MeV`：`301.28034780382563` -> `301.2802511220666`；
  - `rho_hadron`：`1.4970942475813287` -> `1.497012893021389`；
  - `rho_quark`：`1.9198525110836209` -> `1.9197509450122145`；
  - `area_residual`：`4.0838893641649034e-5` -> `3.3141877009250365e-8`。
- CEP 行从旧的 found/finiteness 记录变为当前三态合同的 ambiguous/not-found/NaN 记录。
- crossover 行按当前最大 confirmed boundary chemical potential 重新导出；它们仍不是
  regression 的逐点 crossover oracle，测试只约束行数、端点和有效性边界。

行数、CSV schema、boundary 数量、crossover 数量以及 `T=120 MeV` boundary 数值保持稳定；
没有放宽任何数值容差，也没有修改旧 v2 文件。

## 验证

- `julia --project=. scripts/dev/export_phase_pipeline_regression_baseline.jl --output tests/baselines/phase/baseline_phase_pipeline_v3.csv`：通过。
- 使用独立输出路径重复导出：通过，结果与 v3 字节级一致。
- parent baseline focused regression：`46 passed / 3 failed`，失败仅为 T=125 boundary 的
  `rho_hadron`、`rho_quark`、`area_residual`；该失败作为本次 candidate 的直接诊断证据。
- PR 分支上的 v3 focused regression：`49/49` 通过（`1m36.6s`）。
- 统一 regression runner（`REGRESSION_FILES=phase/test_phase_pipeline_regression.jl`）：`49/49`
  通过（`1m55.7s`）。
- `git diff --check` 在提交前执行并通过。

## 非目标与残余风险

- 不删除或改写 v2，不晋升 phase reference，不启动 transport production。
- 本地 focused regression 通过不替代 `nightly-full-suite`；完整 CI 结果仍是合并 gate。
