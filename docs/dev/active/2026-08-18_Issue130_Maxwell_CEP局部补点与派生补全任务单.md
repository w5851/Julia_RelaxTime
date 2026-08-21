# Issue #130：Maxwell CEP 近端局部补点与派生补全任务单

状态：review；276-target numerical expansion 与同源 solver-free aggregate replay 均已完成并通过
`expansion_candidate`；strict/derived/render 三层已生成，当前仍为 diagnostic-only，等待作者审核，
不晋升 phase-reference。
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

## v7 冻结后的扩展边界

v7 crossover 派生包已在 commit `46f10270` 冻结；它只复制 Maxwell native rows，不能
替代本 route 的真实补点。现有
`.github/workflows/pnjl-issue130-maxwell-cep-local-pilot.yml` 的矩阵和 collector
严格限定 11 个 `pilot_candidate`，不能直接承载 276 个 `input_incomplete` 目标。
因此完整 Maxwell 补点必须新增独立的 versioned expansion workflow/runner/collector，
并保持 pilot artifact 不可变、使用同一 calculation SHA。按 preflight 的 641 个基础
rho key 计，276 个缺失 boundary 目标至少需要 `176,916` 个基础 rho keys；局部
midpoint/refinement 还需按 cap 另计，wall time 尚未测量。扩展 action 必须先通过
solver-free matrix/hash/provenance 检查，并设置 failed-only 重跑和成本止损；在该
versioned workflow 合并前不触发数值扩展。

当前 expansion contract 已落在：

- workflow：`.github/workflows/pnjl-issue130-maxwell-cep-local-expansion.yml`；
- runner：复用 `pnjl_issue130_maxwell_cep_local_pilot.jl`，默认 pilot 仍严格为 11 个
  `pilot_candidate`，expansion 显式使用 `input_incomplete`；
- collector：复用 `collect_pnjl_issue130_maxwell_cep_local_pilot.py`，通过
  `--selection`、`--schema-version`、`--expected-count` 和 `--candidate-verdict`
  隔离两种 artifact；expansion schema 为
  `pnjl_issue130_maxwell_cep_local_expansion_v1`，目标数固定为 276；
- expansion 已完成：numerical run `32354095831`、同源 aggregate replay run
  `32450188629`。旧 replay manifest 曾因 workflow 未持久化 `SOURCE_RUN_ID` 而显示
  `null`；该 provenance 合同已在本分支修复，新的 solver-free replay `32451053476`
  复用同一 source run，完成后以带 `source_run_id` 的 manifest 为准。历史 replay
  不覆盖 numerical evidence。

## 首次 expansion dispatch 诊断（2026-08-20）

- run `32349984548` 使用 merge SHA `64c9b337d9bcbd0a995199e935a336509e227d11`
  触发，prepare 正确验证了 276 个唯一 `input_incomplete` target 和 `176,916`
  个基础 rho keys。
- numerical matrix 没有生成：workflow 将 JSON array 直接传给
  `strategy.matrix`，Actions 要求 object，因而 numerical job 未启动；aggregate
  只得到 0 个 target artifact，`materialized_target_count=0`，verdict 为
  `pilot_inconclusive`。该 run 没有调用 PNJL solver，不能作为数值失败或 partial
  artifact，也不得在 failed-only 恢复中复用。
- repair 边界仅把 prepare 输出改为
  `{ "target_id": [ ... ] }` matrix object，并增加 regression assertion；不改变
  target list、calculation SHA、runner、Maxwell、geometry、容差或 pilot artifact。

## 第二次 expansion dispatch 诊断（2026-08-20）

- run `32352228986` 使用 repair merge SHA
  `f6832544ef1e7bedfd66ddb66829d486a37984cc` 触发；prepare 已成功验证 276 个
  target，但 numerical matrix 仍未展开。
- 根因为 GitHub Actions 单个 matrix 的组合上限为 256，而本 expansion 有 276
  个 target；aggregate 只下载到 0 个 numerical artifact，run 失败。该 run 同样
  没有调用 PNJL solver，不能解释为 solver/Maxwell 失败，也不能作为 partial
  numerical artifact 复用。
- 后续 repair 将同一 target list 按稳定顺序拆为两个各 138 项的 matrix job，
  aggregate 同时等待两个 job；target、artifact 名称、failed-only 语义、calculation
  SHA、runner、Maxwell/geometry 合同和容差均保持不变。

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

## 2026-08-21：Maxwell expansion 与三层派生包收口

- numerical expansion run：`32354095831`；276/276 target 成功，所有最终 rho-mu
  曲线 finite/converged、candidate 唯一且 geometry certificate 通过。冻结证据目录为
  `D:\Desktop\Julia_RelaxTime_issue130_artifacts\maxwell_cep_local_expansion_32354095831_frozen_v1_final`，
  其 `freeze_manifest.json` 保存完整文件 hash inventory；`reference_write=false`、
  `oracle_labels_consumed=false`，数值源的 `solver_called=true`。
- 同源 aggregate replay：旧 run `32450188629` 已 materialize 276/276，但 manifest 的
  `source_run_id` 因 workflow 环境变量未持久化而为 `null`，只作为历史 replay 诊断保留；
  provenance repair 后的 replay run 为 `32451053476`，应验证 `source_run_id=32354095831`、
  `run_mode=aggregate_replay`、`solver_called=false`、无 missing/unexpected/error。
- 三层输出：
  `docs/analysis/pnjl/phase_reference/issue130_phase_reference_layers_v1/`。
  `strict_reference_v1` 保留 C2/v6 unresolved 并追加 276 个真实 Maxwell 行；
  `derived_reference_v1` 使用 `xi=-0.5:0.00625:0.5`，仅在共同 support 内插值，派生行
  标记 `interpolated_noncertified`；`phase_surface_render_v1` 从结构化表生成 no-triangulation
  render 和 coverage mask。三层均保持 `reference_write=false`。
- 当前 verdict：`awaiting_author_review`。作者审核 strict/derived/render 数据、coverage、
  代表图和 provenance 后，才另行执行 phase-reference promotion gate；RS transport 继续暂停。
