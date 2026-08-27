# Issue #130 RS transport phase-reference parity v2

这是 PR #272 合并 SHA `3b19246fb911be4a2efa75fbe14fcb9a793ca256` 上的 solver-free consumer smoke 与
versioned result-import integrity evidence。它不重新运行 equilibrium solver，
也不把 diagnostic-only 数值结果晋升为默认 runtime。

## 固定输入

- calculation SHA: `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- workflow head: `22874505877491754eed27519ad8a7b871c82571`
- source run: `32354095831`; aggregate replay: `32451053476`
- candidate: `data/reference/pnjl/issue130_phase_reference_v1`
- legacy snapshot: `data/reference/pnjl/legacy_phase_reference_v1`
- merge SHA: `3b19246fb911be4a2efa75fbe14fcb9a793ca256`

## Smoke result

`decision.json` verdict: `rs_consumer_smoke_pass_diagnostic_only`。
runtime 使用 strict candidate certified-only view，并保留逐键 legacy fallback；
`legacy` 模式验证显式 `--phase-reference-mode legacy` rollback；diagnostic 模式
只用于全 candidate 行审计。两套 dry-run 的 phase anchor 明确使用
`reference_interpolation`，避免 dry-run 触发 coexistence solver。

## Imported result integrity

`result_import_inventory.csv` 验证 mode-A `910` scan /
`38220` diagnostic rows 和 mode-B
`909` scan / `38178` diagnostic rows，
包括 manifest file hashes、figure source hashes、direct-coexistence side-point
合同、以及旧 `prod_v1` tree hash。两套 result 仍是
`imported_candidate` + `diagnostic_only`，不切换默认 runtime。

## Provenance note

`consumer_source_smoke_runner_raw.json` 保留 runner 的原始 JSON；v2 的
`consumer_source_smoke.json` 仅规范化 legacy snapshot locator，并记录
`legacy_root_normalized=true`。`raw_sha256` 与规范化后的 hash 见 manifest。

## Boundary

本包支持 source-resolution、fallback/rollback 和导入完整性；不代表 RS 数值
parity、全域 production convergence、旧 reference 删除或新的 solver run。
