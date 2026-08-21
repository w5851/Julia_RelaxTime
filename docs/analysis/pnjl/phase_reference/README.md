# Phase-reference 决策证据层

本目录收纳 Issue #130 phase-reference 线的汇总、限定证据和人工 overlay 决策审计。三个包回答的问题不同，保留各自的 `manifest.json`、表格、图和 verdict；本目录只提供逻辑入口，不合并或改写下游 C1/C2/CEP/Maxwell evidence。

## Evidence Packages

| 逻辑角色 | 当前路径 | verdict / 用途 |
| --- | --- | --- |
| current state freeze | [`phase_reference_current_state_freeze_v1/`](phase_reference_current_state_freeze_v1/) | `diagnostic_state_frozen_promotion_blocked`；冻结已有证据状态 |
| limited evidence audit | [`phase_reference_limited_evidence_audit_v1/`](phase_reference_limited_evidence_audit_v1/) | `raw_curve_coverage_complete_diagnostic_only`；审计 unresolved raw curves |
| manual overlay promotion audit | [`phase_reference_manual_overlay_promotion_audit_v1/`](phase_reference_manual_overlay_promotion_audit_v1/) | `blocked_manual_overlay_inconclusive`；记录人工 overlay 不能覆盖自动 gate 的边界 |
| Issue #130 derived layers | [`issue130_phase_reference_layers_v1/`](issue130_phase_reference_layers_v1/) | `awaiting_author_review`；strict Maxwell expansion、uniform-xi derived tables 和 no-triangulation render |

## Reading Order

1. 先读对应包的 `README.md` 和 `manifest.json`，确认生成时间、输入来源和 verdict。
2. 再读 `decision.json`、`execution_log.md` 和 `tables/claim_ledger.csv`，区分 observation、diagnostic candidate 与 promotion gate。
3. 最后查看表格和图 manifest；`full_hybrid_candidate`、raw S-shape 或人工 bracket 都不等于 phase-reference promotion。

## Boundary

- 这些目录只保存诊断/决策证据，不写入 `data/reference/**`，不启动 RS transport production。
- `raw_curve_archive_v1/` 是独立的外部归档指针，不属于本分组的汇总包。
- `issue130_phase_reference_layers_v1/` 仍是诊断/作者审核候选；其中 derived/render 层不能替代 strict 证书，也未写入 `data/reference/**`。
- 本次 namespace 迁移只改变物理目录位置；包内 manifest、checksum、执行日志和生成时路径保持原样，以保留历史 provenance。
