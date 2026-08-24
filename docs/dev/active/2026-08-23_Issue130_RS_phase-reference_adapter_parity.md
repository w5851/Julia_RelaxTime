# Issue #130：RS transport phase-reference adapter parity 与限定重验

状态：active。runtime-switch PR #260 已合并到 `main@1ccf29310fb20c30bcd154f0b4966e25a7565225`。
本任务承接 RS transport 的 consumer parity；它不删除旧 reference，不覆盖旧 transport
production，也不把 solver-free smoke 当作数值 production 通过。

## 固定输入与边界

- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- candidate source run：`32354095831`
- aggregate replay：`32451053476`
- candidate：`data/reference/pnjl/issue130_phase_reference_v1/`
- legacy fallback：`data/reference/pnjl/{boundary.csv,cep.csv,crossover_dense.csv,spinodals.csv}`
- runtime view：strict candidate 的 certified-only 行，缺失/不合格键逐键 legacy fallback
- rollback：`--phase-reference-mode legacy`
- 本机不调用 equilibrium solver；不启动新的 C0/C1/C2、M4 production 或 transport 全量扫描

## 阶段与验收

- [x] 在合并 SHA 上完成 runtime promotion/consumer smoke；candidate/fallback 计数、candidate
  manifest SHA、`solver_called=false` 和 explicit legacy rollback 均有记录。
- [x] 通过实际 `run_gap_transport_scan.jl`/phase-guided reference resolution 入口执行
  `runtime`、`diagnostic`、`legacy` 三模式 source smoke；检查 `mu_q`/`mu_B` 映射、phase kind、
  fallback 语义和 plan anchor。
- [x] 固化 portable evidence：输入哈希、raw dry-run provenance、source summary、parity table、
  claim ledger 和停止边界。
- [ ] 作者审核 solver-free parity evidence。
- [ ] 经作者授权后触发限定 numerical RS pilot；至少比较 candidate runtime 与显式 legacy
  rollback 的 phase anchor、失败/非有限点、transport 输出字段和运行成本。pilot 失败时停止，
  保留 diagnostic-only、fallback、rollback 和旧 production。
- [ ] phase-guided/Paper P1 等消费者完成不再直接依赖旧文件的引用审计后，才创建独立
  `old-reference retirement` PR。

## 当前证据

见：
`docs/analysis/pnjl/phase_reference/issue130_rs_transport_runtime_parity_v1/`。

- runtime view：`certified_candidate_with_legacy_fallback`
- candidate strict：boundary `7162`（certified `3091`）、crossover `1343`、CEP `93`（certified `90`）、
  spinodal `6886`
- fallback：boundary `35`、crossover `315`、CEP `1`、spinodal `31`
- `xi=0, muB=150 MeV, alpha_T=1` 的 phase-guided crossover anchor：
  candidate runtime `198.2066818967471 MeV`，legacy `198.12550509377107 MeV`，差异
  `0.08117680297604579 MeV`。这是 reference 输入差异，不是 RS transport 数值 parity 结论。

## CI 留痕

PR #261 最新提交 `31ab0490a68be02d2f8cde7739a79128b3d3fd23` 的 CI 已于 2026-08-23
全部通过（12 项，无失败或跳过）：unit-smoke、integration-core、Julia/Python script smoke、
PNJL benchmark，以及 task-ledger、active-docs、docs consistency、data-output-path、
legacy-switch、dependency audit 和 advisory governance。

这只证明代码、证据 schema 与 solver-free consumer smoke 合同通过；不改变 numerical pilot
仍需作者审核和单独授权的边界。

## 非目标与回退

- 不删除或重写 `boundary.csv`、`cep.csv`、`crossover_dense.csv`、`spinodals.csv`。
- 不把 `unresolved`、`ambiguous`、`interpolated_noncertified` 行升级为 candidate certified。
- 不移除逐键 fallback，不移除显式 legacy rollback。
- numerical pilot 或后续 consumer 审计失败时，维持旧 reference/旧 production，提交具体失败字段和
  provenance；不放宽物理/数值容差。
