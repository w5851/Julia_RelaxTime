# Issue #130：phase-reference runtime switch 任务单

状态：active。PR #259 已完成 Julia/Python consumer adapter migration；本任务只处理
独立的 runtime source switch，不改 solver、Maxwell、C0/C1/C2 artifact、物理容差或旧
reference 文件。

## Runtime 合同

- 默认 source：`candidate`，layer：`strict`。
- candidate 只暴露 `certified=true` 的行；`unresolved`、`ambiguous`、
  `interpolated_noncertified` 永远保留在 candidate evidence 中，但不能进入 runtime view。
- candidate 缺失的 boundary/crossover/CEP/spinodal key 才由 legacy source 补位，补位行标记
  `source_layer=legacy_fallback`、`status=legacy_fallback`。
- `--phase-reference-mode legacy` 是显式 rollback；`diagnostic` 只允许审计 candidate，
  `render` 不得成为 runtime 输入。
- effective config、run manifest 和 phase-guided sidecar 记录 source、runtime view、fallback
  reason、candidate manifest SHA 和逐表 candidate/fallback 计数。

## 验收

- Julia fixtures：candidate certified/unresolved、逐键 fallback、legacy rollback、mu_q/mu_B
  parity、NaN/Inf/duplicate key 和默认 CLI 解析。
- Python fixture：与 Julia 相同的 certified-only/fallback schema 语义。
- solver-free：`run_gap_transport_scan.jl` 和 phase-guided CLI 的 default/legacy/diagnostic
  dry-run 只检查 source resolution，不调用 equilibrium solver。
- focused unit/integration、task-ledger、docs consistency、`git diff --check` 全部通过后，
  创建独立 runtime-switch PR；PR 合并并经作者审核前不称 candidate 已正式晋升，也不启动
  RS transport 或旧 reference retirement。

## Provenance

- Base merge：PR #259，`960dc751f5efaface38b0d56db56ba577f67fd2d`。
- Candidate root：`data/reference/pnjl/issue130_phase_reference_v1/`。
- Legacy rollback files：`boundary.csv`、`cep.csv`、`crossover_dense.csv`、`spinodals.csv`。
