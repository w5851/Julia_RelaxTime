# `docs/analysis` Metadata Repair v1

状态：`accepted`

日期：2026-08-20

## 范围与判定

本批只修复以下四个诊断包的历史输出 metadata：

- `pnjl/phase_reference/phase_reference_current_state_freeze_v1`
- `pnjl/phase_reference/phase_reference_limited_evidence_audit_v1`
- `pnjl/phase_reference/phase_reference_manual_overlay_promotion_audit_v1`
- `relaxtime/phase_guided_transport/phase_guided_transport_p128_xi001_analysis`

修复前的全部 40 个 mismatch 均判定为 `stale_metadata`：

| package | manifest mismatch | checksum mismatch | classification |
| --- | ---: | ---: | --- |
| `phase_reference_current_state_freeze_v1` | 8 | 9 | `stale_metadata` |
| `phase_reference_limited_evidence_audit_v1` | 11 | 0 | `stale_metadata` |
| `phase_reference_manual_overlay_promotion_audit_v1` | 10 | 0 | `stale_metadata` |
| `phase_guided_transport_p128_xi001_analysis` | 2 | 0 | `stale_metadata` |

历史基线为 phase-reference 的首次入库提交 `a564d3288ba35913172894c3ccc47db6f90d68d7` 和 transport v1 的首次入库提交 `5f75efab086bccac95eb07e84b7c5c15c8f011ba`。当前文件与各自基线逐字节一致；没有发现 `newline_or_serialization_drift`、`content_change` 或需要 `author_check` 的 mismatch。

## 修改内容

- 更新三个 phase-reference `manifest.json` 的 `output_files` SHA-256。
- 更新 current-state 包的 `checksums.sha256`，包括修复后的 `manifest.json` hash。
- 更新 transport v1 root manifest 的 README/plot manifest hash 与 bytes，并将 21 个 output 路径规范化到当前 namespace。
- 将 transport v1 `figures/plot_manifest.json` 中的 19 个 figure/input-table 路径规范化到当前 namespace；图像 hash 和 bytes 未变。

未修改 CSV、PNG、README 科学内容、verdict、source provenance、输入文件 hash、claim ledger、solver、Maxwell、reference 或 transport 数值产物。phase-guided transport v2 的 52 个 output 保持范围外且未修改。

## 审计证据

- `pre_repair_snapshot.json`：保存四个包的文件 inventory、当前字节 SHA-256、旧 manifest/checksum 原文、历史基线和逐项分类。
- `post_repair_audit.json`：修复后严格审计；四个包均为 `0 mismatch / 0 missing / 0 extra`。
- pre-repair snapshot SHA-256：`b4497989a194e4b36d392f4dc38036e2fddd38f7c7e5aad15bb9ef5d0dadcd31`
- post-repair audit SHA-256：`995e5c33938cacb577b62792e5e806afa66be1edb34338646a91455e132936e0`

独立不变性检查结果：四个包的非 metadata payload `payload_invariant=True`；去除输出 hash registry 后 manifest 语义 `semantic_manifests_invariant=True`；四个 claim ledger 文件 `claim_ledgers_invariant=True`；transport figure manifest 除路径 namespace 外内容等价。

## 验证

- `python scripts/dev/audit_analysis_metadata.py --strict`：通过。
- 四个目标包的 JSON 解析：通过。
- `git diff --check`：通过。
- `julia --project=. scripts/dev/check_task_ledger.jl`：通过。
- `julia --project=. scripts/dev/check_task_ledger.jl --preflight --track analysis-docs-cleanup`：通过；dirty paths 均属于本批 allowlist。
- `julia --project=. scripts/dev/check_docs_consistency.jl`：通过。
- `julia --project=. scripts/dev/check_active_docs_governance.jl`：通过。
- `julia --project=. scripts/dev/check_script_entrypoints.jl`：通过，检查 203 个 Julia 文件。
