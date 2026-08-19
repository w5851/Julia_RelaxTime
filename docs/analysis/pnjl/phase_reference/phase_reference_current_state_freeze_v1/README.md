# Issue #130 当前结果冻结 v1

生成时间：`2026-08-18T06:19:13Z`。本包是当前 phase-reference 讨论的状态快照，目的
是把已有 C0/C1/C2、full-shadow、C2 grid、raw-curve audit 和人工 CEP overlay
放在同一个可追溯入口下。它不重跑 solver，不修改历史 evidence，不写入
`data/reference/**`，也不把诊断候选升级为生产 reference。

## 固定身份

- 当前工作树：`D:\Temp\julia_relaxtime_phase_reference_limited_evidence_audit`
- 当前分支：`codex/issue130-phase-reference-limited-evidence-audit`
- 当前 HEAD：`e25047e871cca730e14a9783aee63fbff42a76dc`；工作树在生成时为 dirty，原有改动保持不动。
- 当前数值 calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`。
- 最新 Stage-C C2 source run：`31862752226`；postprocess SHA：`fd359e792a89beb5ab12349bba761dc58ee16761`。
- full endpoint-local shadow aggregate：`31714535418`；结果为 `full_hybrid_candidate`。
- 历史 C0/C1/C2 comparator 使用 SHA：`ffa816df0a145f73d7490db1ed9ff10c92e017a4`，仅作历史对照。

## 已固化事实

1. 最新 C2 primary rows：boundary `6886`、spinodals `6886`、crossover `1488`、
   CEP `93`，输入文件均为 finite，主表 `converged=false` 为零且声明的主键无重复。
   C2 phase diagnostics 的状态计数为 confirmed first-order `6886`、confirmed
   monotone `2636`、ambiguous `359`。
2. `phase_grid_convergence` 共 `22791` 行，其中 `5424` 行 unresolved。它是逐层
   refinement 事件表，不是 Maxwell boundary 表；同一 `(xi,T)` 出现多层记录是设计语义，
   因此 unresolved 行只能保留为诊断，不能直接绘制成有效 Maxwell 面。
3. 独立 raw-curve audit 覆盖 `763` 个 unresolved 行、去重后 `761` 个坐标；每个坐标
   有 `1281` 个 rho 点且 finite/converged。raw 曲线中的 `+→−→+` 只是拓扑观察，
   不等于 Maxwell 面积或 geometry certificate。
4. 九个 targeted regression 点的自动 overlay 仍是 `hybrid ambiguous` 对
   `oracle confirmed_first_order`，`classification_match=false`；作者对临时曲线的
   视觉判断只记录为 `accepted_diagnostic_only`，没有覆盖自动 gate。
5. 三个固定 T 方向 CEP 切片已经人工审核：`xi=0.125`、`0.39375`、`0.5` 的
   bracket 均为 `0.0625 MeV`，状态为 `closed_by_author_review`。这关闭的是固定切片
   的诊断 bracket，不是单值 CEP 写回。
6. μ 方向 CEP endpoint refinement 尚未运行；下一步仍是 solver-free target-list/cost
   preflight。不得把已被排除的 `mu_q > mu_CEP` response peak 当作 crossover 锚点。

## 门禁结论

`full_hybrid_candidate` 说明 endpoint-local shadow 的自动候选 gate 通过，不能等同
于 phase-reference promotion。当前 strict promotion 仍被 grid unresolved 和九点
自动分类不一致阻塞；人工 overlay 只增加诊断事实，不改变 C2 artifact。RS transport
继续依赖 phase-reference promotion，尚未启动。

详见：

- `tables/c2_summary.csv`
- `tables/gate_status_matrix.csv`
- `tables/input_validation.csv`
- `tables/cep_manual_overlay.csv`
- `tables/targeted_regression_status.csv`
- `tables/claim_ledger.csv`
- `figures/plot_manifest.json`
- `checksums.sha256`（生成后写入）

## 明确不在本冻结包中的内容

- 不重跑 C0/C1/C2，不启动 μ numerical pilot、reference write、C3/O1 或 transport。
- 不覆盖 C2 的 unresolved、CEP bracket、Maxwell candidate 或已有图像。
- 不把插值/视觉闭合图当作 strict reference；后续 derived completion 必须单独版本化并
  标明 `interpolated_noncertified`。
