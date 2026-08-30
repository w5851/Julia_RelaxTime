# Issue #130 RS transport phase-reference parity v1

这是合并 SHA `1ccf29310fb20c30bcd154f0b4966e25a7565225` 上的 solver-free consumer/parity evidence。

- calculation SHA: `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- source run: `32354095831`; aggregate replay: `32451053476`
- runtime view: `certified_candidate_with_legacy_fallback`
- runtime 首选 strict candidate 的 certified-only 行；缺失或不合格键逐键由 legacy fallback 补位。
- `legacy` 是显式 rollback；`diagnostic` 仅用于审计，不是 runtime 输入。
- 本包只运行 source-resolution 和 phase-guided plan dry-run；`solver_called=false`，没有输运系数。

## 结果

candidate 与 legacy 的相变锚点数值不同，但字段/单位/phase-kind 语义一致；差异是 reference 输入差异，不能被称为 RS 数值收敛。
下一步是作者审核本包后，单独触发限定 numerical RS pilot；在该 pilot 通过前保留旧 reference、fallback、rollback，不创建 retirement PR。
