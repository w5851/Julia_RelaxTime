# Issue #130 RS old-reference retirement audit v1

这是 RS transport 旧 `prod_v1` 的 solver-free retirement audit。审计只读取仓库内已入库的 result、figure、registry、manifest、代码和文档；不调用 Julia/equilibrium/transport solver，不改写或删除任何数值产物。

- audit schema: `issue130_rs_old_reference_retirement_audit_v1`
- audit repo HEAD: `ac0674e249176a3453257105bab3b2382e409581`
- calculation SHA: `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- workflow head: `22874505877491754eed27519ad8a7b871c82571`
- verdict: `retirement_audit_pass_retain_legacy`
- solver_called: `False`
- production_write: `False`

## 范围和结论

本审计的 old reference 指 RS phase-guided transport 的 `first_canonical_v2...prod_v1` mode-A/mode-B。更早的 `first_canonical_v1`/`xi005` 树另列为历史证据，不自动纳入当前 runtime fallback 的删除范围。PNJL phase-reference 的 `legacy_phase_reference_v1` 是另一条已经完成 canonical-path retirement 的链，本包不修改它。

当前树族计数：`{'current_prod_v2': 4, 'historical_prod_v1': 10, 'legacy_prod_v1': 4}`。旧 `prod_v1` 结果的 fallback tree hash 与两套当前 `prod_v2` manifest 的记录逐模式核对；详见 `fallback_smoke.csv`。
较早历史树的 manifest hash mismatch 计数为 `26`；这些只作为历史 provenance warning，不阻断当前 fallback tree 的完整性门禁，详见 `reference_inventory.csv`。

审计建议是：旧 `prod_v1` 可以继续评估从 active/canonical 默认路径退出，但必须保留不可变、带哈希的 legacy snapshot 与显式 rollback；本包不授权物理删除。

## Consumer smoke

默认/回退静态合同检查：`5/5` 通过；旧 reference 的 active consumer occurrence 为 `0`。这不是数值 parity 或 numerical convergence 证明。

- 默认分析入口应解析 current `prod_v2`。
- legacy 只通过显式 versioned snapshot/fallback 访问。
- `production_registry.json` 的 approved/superseded 状态与 raw manifest 内的导入/runtime 状态分开记录；`metadata_review_required` 不能直接当作数据损坏。

## Retention boundary

| 层 | 处理 | 删除状态 |
| --- | --- | --- |
| current `prod_v2` | 保留 canonical | 不适用 |
| current fallback `prod_v1` | 可退出 active/canonical，迁入 versioned legacy snapshot | 未授权，不删除 |
| older `prod_v1`/xi005 | 保留历史和收敛证据 | 需独立历史依赖审核 |

## Evidence files

- `reference_inventory.csv`: 每棵 result/figure tree 的文件数、字节数、tree/manifest hash、schema/finite/duplicate 状态。
- `registry_consistency.csv`: registry 与各 tree manifest 的状态配对和 superseded 关系。
- `consumer_reference_map.csv`: tracked text 中旧 token 的逐行引用、分类和 review role（active consumer / fallback-or-rollback / historical 等）。
- `fallback_smoke.csv`: current `prod_v2` 到旧 `prod_v1` 的逐模式 hash/fallback/rollback 检查。
- `retention_deletion_allowlist.csv`: 只读的保留、路径退役和物理删除边界。
- `rollback_plan.md`: 后续 implementation PR 的回退顺序；没有删除命令。
- `claim_ledger.csv`: 可支持结论与未声明结论。

## Stop conditions

任一 active old consumer、current/fallback tree hash 不匹配、未能解析 versioned snapshot、或历史证据依赖未厘清时，停止 implementation PR。较早历史树的已知 manifest drift 必须原样记录，不能改写成当前 fallback 通过。不得删除旧树、不得重跑 solver，也不得把本包写成 RS 数值收敛证据。
