# Issue #130 PNJL accepted-primary / legacy retirement audit v4

这是 solver-free 的合同审计：通过 production-parity Python adapter 读取 v2 `accepted` 与显式 `strict`，
再按真实 semantic key 对照 byte-preserving legacy snapshot。它不调用 PNJL solver、不写 reference、不删除文件。
repo HEAD：`7b02d84d3c38c60500e82b769eb0f9c22de39c9b`；calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`。

## Verdict

- `legacy_physical_deletion_candidate`
- accepted primary runtime ready：`True`
- runtime legacy fallback rows：`0`；rollback enabled：`False`
- active legacy path contracts pending migration：`0`
- physical deletion eligible：`True`
- next action：prepare a separate exact allowlist deletion PR and request author review

## 三个不混淆的覆盖问题

1. accepted 是否可作为 runtime primary：由 package promotion、adapter eligibility 和 strict 显式路径验证。
2. legacy 是否仍会被 runtime fallback/rollback 调用：由 adapter API 和 runtime mode 扫描验证；当前合同计数应为 0。
3. 物理删除前是否仍有旧路径文本：由 consumer matrix 逐项迁移；这与数值 exact-key 覆盖是不同门禁。

## Semantic-key coverage（仅诊断）

accepted 与 legacy 的 exact key 不要求相同。`neighbor_coverage.csv` 记录同 ξ 最近 accepted 支持，
但 nearest support 不能伪装为 exact legacy key；consumer 应使用 accepted 自己的插值/最近 ξ 语义。

| table | accepted | strict certified | legacy | exact accepted | physical missing | above-CEP excluded |
|---|---:|---:|---:|---:|---:|---:|
| boundary | 12537 | 7162 | 48 | 14 | 34 | 0 |
| crossover | 3135 | 1343 | 336 | 21 | 190 | 125 |
| cep | 161 | 90 | 11 | 11 | 0 | 0 |
| spinodals | 11989 | 6886 | 57 | 26 | 31 | 0 |

## Path-retirement boundary

legacy snapshot 仍保留，直到 active path contracts 迁移完成并由独立 physical-deletion PR 审核。
历史 audit、snapshot manifest 和 contract tests 可以保留；它们不等于 runtime fallback。

## Provenance

- accepted package manifest SHA：`3e37d6812c33c4dd15469391247768fe536dc3185fd54543cdaba2713c520be1`
- accepted layer manifest SHA：`e053824282c556f252e952163a534c35ddad554452ade6130c2411e11566f60d`
- accepted package tree SHA：`a0a76331e4a1409f7b89deeffd76bb0eb77097068c846b39af125678912e0c30`
- legacy retirement manifest SHA：`1bb479e1876ac94157d5f2d04e435ff60fcb7e4e2615f169c22680d7aec1ce1d`
- legacy snapshot tree SHA：`c128ef6358a5813533fc5a9726047585a750a4421f162293564bea8e363764e6`
- solver_called：`false`；reference_write：`false`；runtime_consumption：`false`（本审计自身）。
