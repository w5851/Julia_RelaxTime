# Issue #130：Maxwell CEP 近端局部补点与派生补全任务单

状态：active；11-target numerical pilot 与同源 aggregate replay 均已完成并通过
`pilot_candidate`；当前作为 v7 冻结后的下一步 Maxwell 输入 route，证据仍为
diagnostic-only，不晋升 phase-reference。
父任务为 `issue130-phase`。这是与 crossover μ endpoint refinement 分开的
required follow-up，只处理 Maxwell 侧在 CEP 附近的 support/geometry 缺口。

## 目标

在固定 C2 v5 证据上，先用 solver-free 方式建立 Maxwell 近 CEP 的 target list
和成本 envelope；只有预检通过并获得单独授权后，才对少量代表性 `xi` 运行
numerical pilot。目标是判断缺失来自 rho geometry/refinement 还是端点采样，
不是把 unresolved 记录直接改成通过。

## 固定输入与预检结果

- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`。
- v5 证据：`docs/analysis/pnjl/c2_phase_surfaces_diagnostic_v5_no_triangulation/`。
- 目标表：`maxwell_surface_point_status.csv`、`grid_unresolved_diagnostics.csv`、
  `crossover_maxwell_endpoint_separation.csv`。
- preflight 输出：`docs/analysis/pnjl/issue130_endpoint_refinement_preflight_v1/maxwell_local/`。
  窗口为 `T <= T_CEP_bracket_low` 且距离不超过 `8 MeV`；共 760 个候选，484 个
  已有 boundary 行，276 个为 `input_incomplete`，代表性 xi 共 11 个 pilot 目标。

## 预检与 numerical pilot 合同

- 预检按固定 `(xi,T)` 识别 CEP 下方 unresolved Maxwell 行；缺少 boundary 行的
  目标只保留为 `input_incomplete`，不得自动进入 numerical。
- 每个数值目标必须重新取得完整 `rho-mu` 曲线，执行公共 Maxwell、三态和 geometry；
  candidate 必须唯一，三个交点 finite/有序，且不放宽 position/density/area 门禁。
- 成本只报告 `rho=0:0.00625:4` 的 641 点单位；wall time 未测量。预检通过不等于
  Maxwell certificate 或 production feasibility。pilot 额外记录最多 12 个当前
  crossing bracket midpoint，
  并按每级 strict candidate 结果计费。

## Numerical pilot implementation

- workflow：`.github/workflows/pnjl-issue130-maxwell-cep-local-pilot.yml`；
  runner：`scripts/analysis/pnjl_issue130_maxwell_cep_local_pilot.jl`；
  collector：`scripts/analysis/collect_pnjl_issue130_maxwell_cep_local_pilot.py`。
- 固定 11 个 `pilot_candidate` target；每个目标先完整计算
  `rho=0:0.00625:4`，随后最多补 12 个当前三交点 bracket midpoint；不读取 oracle
  标签，不写 phase-reference。
- focused CI 已通过并合并；使用 calculation SHA
  `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48` 的 numerical pilot run 为
  `32222254605`，workflow/postprocess SHA 为
  `51c93cf0111415b35bb199376c358782c0f5a2f4`。11/11 target 均生成 653 点完整
  曲线、唯一三交点和 geometry convergence；总 unique solves 为 7183，targeted
  additions 为 132，solver failure/nonfinite/retry/fallback 均为 0。数值 artifact
  已下载到 `D:\Desktop\Julia_RelaxTime_issue130_artifacts\maxwell_cep_local_pilot_32222254605`。
- numerical run 后的初始 aggregate verdict 为 `pilot_inconclusive`，历史原因是
  `target_summary.json` 缺少 calculation/workflow SHA；对应 `provenance.json` 与
  `manifest.json` 的 SHA 正确。该历史失败已由后续 provenance repair 收口；summary
  缺失字段只在两个身份源一致且符合输入时允许 fallback，summary 存在但错误、或
  两个身份源不一致时仍失败。
- PR #239 merge SHA 为 `ebf8e9199ff81bbb31df7f38cf0368a0040c103a`；其首次 replay
  run `32225731719` 的 `pilot_inconclusive` 仅由 source/postprocess SHA 角色混用造成，
  不涉及 solver、曲线、Maxwell 或 geometry。
- PR #240 merge SHA 为 `869e1f0f1ae5bd490cacfb6066fd27205f054100`。正式 replay
  run `32237794907` 复用数值源 run `32222254605`，以 source workflow SHA
  `51c93cf0111415b35bb199376c358782c0f5a2f4` 校验 numerical artifacts，并单独记录
  `postprocess_sha=869e1f0f...`；`solver_called=false`、`reference_write=false`、
  `oracle_labels_consumed=false`。11/11 target materialized，verdict=`pilot_candidate`，
  aggregate errors=0。该 replay 不重跑 solver，Maxwell 证据仍保持 diagnostic-only。

## Derived completion 边界

`certified_layer` 保留原始 C2 证书；`completed_layer` 只能在同一 Maxwell 物理区、
相邻已证实行之间做显式 `interpolated_noncertified` 派生，不能伪造 candidate、
`maxwell_area` 或 strict geometry certificate，也不能自动晋升 phase-reference。

本 companion route 是 Issue #130 三层验收目标中 `strict_reference_v1` 的 Maxwell
输入来源，也是后续 `derived_reference_v1` 和 `phase_surface_render_v1` 的前置依赖。
Maxwell 真实补点完成前，不得用 crossover 派生层的视觉连续性填充 Maxwell；最终
render 只能由统一 ξ 的 derived reference 生成，且必须保留 Maxwell unresolved
mask 和 cell coverage。

## 与 crossover 全 ξ expansion 的边界

Maxwell pilot/replay 已满足其自身的 11-target diagnostic gate，但不授权自动补算
276 个 `input_incomplete` 行，也不改变 crossover μ endpoint refinement 的目标表。
作者已另行授权 crossover 全 ξ expansion；该任务使用相同 calculation SHA、独立
versioned workflow 和独立 aggregate，两个 route 的 target、成本和 verdict 分开记录。

不修改 equilibrium solver、Maxwell、三态规则、endpoint policy 或容差；不重跑 C0/C1/C2。

## v6 派生相图（2026-08-20）

- v6 不直接读取 C2 原始结果重建 surface，而是以 `c2_phase_surfaces_diagnostic_v5_no_triangulation` 的后处理表为唯一 baseline。
- Maxwell boundary、spinodal、CEP bracket、v5 crossover 物理筛选、unresolved 状态和 no-triangulation/native-gap 规则均原样保留。
- 仅叠加 crossover endpoint expansion numerical run `32240898122` 的 solver-free aggregate replay `32255786553`；固定 calculation SHA 为 `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`。
- v6 物化 186 个 endpoint candidate，覆盖 93 个非均匀 ξ 切片；所有点 finite/converged 且 `mu_q <=` 同一 v5 CEP proxy，但仍只表示 diagnostic overlay，不等于 CEP 已闭合或 phase-reference 已晋升。
- v6 产物：`docs/analysis/pnjl/c2_surface_views/c2_phase_surfaces_diagnostic_v6_crossover_overlay/`。Maxwell 276 个 `input_incomplete` 候选仍未补算，本次不扩大 Maxwell 数值范围。
