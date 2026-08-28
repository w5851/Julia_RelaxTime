# Issue #130：RS old-reference retirement audit

状态：accepted；审计 PR #276 已合并到 `main@ad47ad1482574faa4764af3913b8cfa3b6bae53e`，
作者审核已通过，当前进入独立的 path-retirement implementation PR。仍不物理删除旧
`prod_v1`，不调用 solver，不触发 numerical Actions。

## 范围锁定

- 目标是 RS transport 的旧
  `first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1`，覆盖
  mode-A `fixed_muB_phase_scaled` 和 mode-B `fixed_T_sparse_muB` 的 result/figure tree。
- 更早的 `first_canonical_v1`、`xi005` 树只作为历史收敛/出版证据单独保留，不在本审计中自动删除。
- PNJL phase-reference 的 `legacy_phase_reference_v1` 已由 PR #264 完成 canonical-path retirement，
  不属于本任务。
- current `prod_v2` calculation SHA 固定为
  `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`；workflow head 固定为
  `22874505877491754eed27519ad8a7b871c82571`。

## Audit PR

分支：`codex/issue130-rs-old-reference-retirement-audit`

审计脚本：
`scripts/analysis/relaxtime/audit_issue130_rs_old_reference_retirement.py`

审计包：
`docs/analysis/relaxtime/issue130_rs_old_reference_retirement_audit_v1/`

### 已完成

- [x] 逐模式盘点 current `prod_v2`、legacy fallback `prod_v1`、历史 `prod_v1` result/figure tree，记录文件数、字节数、tree SHA-256、manifest SHA-256、scan/diagnostic 行数、重复键、非有限值和失败状态。
- [x] 核对两套 current `prod_v2` manifest 记录的 `legacy_prod_v1_tree_hash` 与实际旧 tree；确认显式 `--phase-reference-mode legacy` 和 fallback 标志仍存在。
- [x] 对 tracked code/config/docs/workflow/manifest 做旧 RS case 引用扫描；区分 active consumer、治理/审计工具、生成产物和历史证据。无 active old consumer 时才允许进入 implementation review。
- [x] 检查默认分析入口使用 current `prod_v2`，并检查 versioned legacy snapshot/source switch 的 solver-free 静态合同。
- [x] 生成 retention/deletion allowlist 与 rollback plan；`old_reference_deleted=false`、`solver_called=false`、`production_write=false`。

### 当前 audit 结果

- verdict：`retirement_audit_pass_retain_legacy`。
- active old reference occurrences：`0`。
- 两个 fallback hash smoke：`pass`；默认/回退静态 consumer smoke：`5/5 pass`。
- 旧 v2 `prod_v1` 可以进入 active/canonical path retirement 的 implementation review，但必须保留带哈希的 versioned legacy snapshot 和显式 rollback。
- 历史 `v1/xi005` tree 有已知 manifest drift（详见 `reference_inventory.csv`），只作为历史 provenance warning，不能改写或冒充当前 fallback 完整性。

## 放行门禁

- [x] 作者审核 `reference_inventory.csv`、`registry_consistency.csv`、`consumer_reference_map.csv` 和
  `retention_deletion_allowlist.csv`。
- [x] 另立 implementation PR，仅迁移旧 v2 `prod_v1` 的 canonical/path 指针、registry 和文档，
  不改 raw bytes、不重算数值；具体执行任务见
  `2026-08-28_Issue130_RS_old_reference_path_retirement.md`。
- [ ] implementation PR 通过 focused CI，并重新执行 default/legacy/rollback solver-free smoke。
- [ ] implementation PR 合并后，旧 raw/figure 只保留 versioned snapshot；任务台账记录新的 mainline SHA。
- [ ] 物理删除必须另立 PR，并在历史依赖、论文复现、rollback 和 snapshot hash 均确认后获得单独授权；
  当前不自动删除。

## 停止条件

任一 active old consumer、current/fallback tree hash 不匹配、versioned snapshot 无法解析、registry
路径语义不一致或历史证据依赖未厘清时，停止 implementation PR；不重跑 solver，不放宽容差，不删除
旧 `prod_v1`。

## Audit merge record

- audit PR：#276，merge SHA `ad47ad1482574faa4764af3913b8cfa3b6bae53e`
- verdict：`retirement_audit_pass_retain_legacy`
- solver/free boundary：`solver_called=false`、`production_write=false`、`old_reference_deleted=false`
- handoff：path-retirement implementation 分支
  `codex/issue130-rs-old-reference-path-retirement`
