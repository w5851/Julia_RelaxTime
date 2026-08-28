# Issue #130：RS old `prod_v1` physical deletion proposal

状态：active；PR 提案分支已应用精确 allowlist 删除，尚未合并。path-retirement
PR #277 已以 `main@74b53b47ebcca2b292cee72f70a70a84b0d2eea5` 合并；本任务只讨论
其后的工作树物理删除，不回写或重跑数值。

## 范围

只删除以下 versioned snapshot：

```text
data/outputs/results/relaxtime/transport/phase_guided/legacy_prod_v1_snapshot_v1/
data/outputs/figures/relaxtime/transport/phase_guided/legacy_prod_v1_snapshot_v1/
```

实际 allowlist 见
`docs/analysis/relaxtime/issue130_rs_old_reference_physical_deletion_v1/deletion_allowlist.csv`。
它包含两个 result 根元数据文件、两棵 result case tree 和两棵 figure case tree，
合计 112 个文件、41,673,063 bytes。更早的 `first_canonical_v1`/`xi005` 历史证据、
当前 `prod_v2`、PNJL `legacy_phase_reference_v1` 均明确排除。

## 物理与回退边界

- 删除前的 tree/manifest hash、文件数和字节数由 deletion manifest 固化；
- Git 恢复引用为 `74b53b47ebcca2b292cee72f70a70a84b0d2eea5`，步骤见
  `restore_commands.md`；仓库没有 RS 专用外部 archive；
- 合并后 RS `prod_v1` fallback/rollback 不再可用，不能继续在文档中声称
  `explicit_snapshot_only`；PNJL phase-reference 的 legacy fallback 不受影响；
- 不调用 solver、不改当前 approved `prod_v2`、不触发 numerical Actions。

## 验收门禁

- [x] PR #277 CI 全绿并已合并；path-retirement 远端/本地分支已清理；
- [x] 仅删除 allowlist 中的 snapshot paths；
- [x] deletion manifest 记录 scope、hash、字节数、Git 恢复引用和排除项；
- [ ] deletion validator、current `prod_v2` integrity 和 focused CI 全绿；
- [ ] 作者审核确认 Git 历史恢复边界和“删除后无 RS fallback/rollback”语义；
- [ ] 物理删除 PR 单独获得合并授权。

## 停止条件

allowlist 以外出现删除、current `prod_v2` 或 PNJL legacy 文件变化、恢复引用不可达、
registry 仍把旧 snapshot 标成可用、或任何 numerical/solver 变化时，立即停止。
