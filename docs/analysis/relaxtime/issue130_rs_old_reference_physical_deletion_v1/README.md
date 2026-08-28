# Issue #130 RS old `prod_v1` physical-deletion proposal v1

本包记录 RS phase-guided transport 旧 `prod_v1` snapshot 的物理删除提案。
它只在物理删除分支中移除工作树副本；当前 PR 尚未合并，`main` 的工作树状态由
path-retirement merge `74b53b47ebcca2b292cee72f70a70a84b0d2eea5` 保留。

## Scope

只删除以下 versioned snapshot：

- `data/outputs/results/relaxtime/transport/phase_guided/legacy_prod_v1_snapshot_v1/`
- `data/outputs/figures/relaxtime/transport/phase_guided/legacy_prod_v1_snapshot_v1/`

这两个根目录下包含两套 mode-A/mode-B `first_canonical_v2...prod_v1` 树和两个
result 根元数据文件。删除前的树 hash、manifest hash、文件数和字节数见
`deletion_allowlist.csv` 与 `deletion_manifest.json`。

明确不删除：

- 当前 approved `prod_v2` result/figure；
- 更早的 `first_canonical_v1`、`xi005` 和 convergence 历史证据；
- `data/reference/pnjl/legacy_phase_reference_v1`（这是 PNJL phase-reference
  的另一条 fallback 链）；
- 任意 numerical solver、Actions artifact 或 production 数值。

## Recovery boundary

仓库中没有为 RS `prod_v1` snapshot 建立独立 Zenodo archive。本提案删除工作树副本，
但保留 path-retirement merge commit `74b53b47ebcca2b292cee72f70a70a84b0d2eea5` 作为
Git 历史恢复引用。恢复步骤见 `restore_commands.md`；合并后若仓库历史被重写或清理，
该恢复边界也会失效，因此这不是新的数值备份。

## Status and gates

- `proposal_status`: `pending_author_review`
- 本分支已应用 `git rm`，但不代表 `main` 已删除；合并仍需单独授权。
- `fallback_available_after_merge=false`、`rollback_available_after_merge=false`：
  这是有意的破坏性边界，不能与旧 snapshot 删除混写成“仍可回退”。
- `solver_called=false`、`production_write=false`：本 PR 只做路径/registry/文档/测试治理。
- 合并前必须通过 deletion validator、当前 `prod_v2` 完整性检查、全套 focused CI，
  并由作者确认 Git 历史恢复边界足够。
