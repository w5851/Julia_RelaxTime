# PNJL Analysis Index

本目录是 PNJL，尤其是 Issue #130 phase-reference 相关的诊断证据集合。各 case 的 `manifest.json`、`decision.json` 和 `AUDIT.md` 是本 case 的边界；本索引只建立研究线和时间线，不把 diagnostic candidate 晋升为 production/reference。

## Issue #130 Logical Groups

### C1 surface views

| 逻辑角色 | 当前路径 | 说明 |
| --- | --- | --- |
| `mu_xi_T` view | [`c1_surface_views/pnjl_c1_mu_xi_T_phase_surfaces_diagnostic_v2/`](c1_surface_views/pnjl_c1_mu_xi_T_phase_surfaces_diagnostic_v2/) | C1 诊断相面，包含 CEP temperature bracket view |
| `xi_T_mu` view | [`c1_surface_views/pnjl_c1_xi_t_mu_phase_surfaces_diagnostic_v1/`](c1_surface_views/pnjl_c1_xi_t_mu_phase_surfaces_diagnostic_v1/) | 同一 source run 的替代坐标投影 |

两者共享 C1 source run 和部分 summary/claim 表，但坐标语义和图 manifest 不同；当前保留为两个独立 view 以维护各自的 provenance。

### C2 phase surfaces and audits

| 逻辑子线 | 当前路径 | 状态/用途 |
| --- | --- | --- |
| surface versions | `c2_surface_views/c2_phase_surfaces_diagnostic_v1/` 到 `c2_surface_views/c2_phase_surfaces_diagnostic_v5_no_triangulation/` | 同一 C2 surface 线的版本化诊断；v5 是当前 native-support/no-triangulation 视图 |
| visual variant | `c2_surface_views/c2_phase_surfaces_diagnostic_v4_visual_closed/`、`c2_surface_views/c2_phase_surfaces_diagnostic_v4_visual_closed_display16/` | v4 的视觉闭合和 display16 变体；不是 reference promotion |
| convergence/blocking | `c2_convergence_audit_v1/`、`c2_blocking_audit_v2/` | C0/C1/C2 比较、阻塞原因和 unresolved 状态审计 |
| review/follow-up | `c2_targeted_manual_review_v1/`、`c2_limited_feasibility_v1/`、`c2_cep_xi05_high_side_extension_v1/`、`c2_manual_bisection_v1/` | 人工复核、输入合同和后续补点，不是同一 surface case 的重复结果 |

### CEP and Maxwell evidence

| 逻辑阶段 | 当前路径 |
| --- | --- |
| local CEP pilot | `cep_narrow_pilot_v1/`、`cep_narrow_pilot_v2/` |
| Stage-C feasibility/replay | `cep_hybrid_stagec_offline_feasibility_v1/`、`cep_hybrid_stagec_certificate_feasibility_v2/`、`cep_hybrid_stagec_extrema_guard_feasibility_v1/`、`cep_hybrid_stagec_tolerance_replay_v1/` |
| cascade/shadow | `cep_cascade_production_shadow_v1/`、`cep_maxwell_endpoint_local_production_shadow_v4_full_31352083775/` |
| endpoint-local contract/shadow | `cep_maxwell/endpoint_local/pnjl_cep_endpoint_local_contract_feasibility_v2/`、`cep_maxwell/endpoint_local/pnjl_cep_endpoint_local_production_shadow_v4/`、`cep_maxwell/endpoint_local/pnjl_cep_endpoint_local_production_shadow_v4_20260813/` |
| Maxwell contracts | `cep_maxwell/maxwell_contracts/pnjl_maxwell_endpoint_candidate_feasibility_v1/`、`cep_maxwell/maxwell_contracts/pnjl_maxwell_endpoint_limit_contract_v1/`、`cep_maxwell/maxwell_contracts/pnjl_maxwell_tolerance_contract_feasibility_v1/` |
| algorithmic feasibility | `criticality_feasibility_v1/` |

这些目录应按阶段串联，而不是把不同 run、不同 contract 和不同 verdict 的 CSV 合并到一个大目录。

### Phase-reference decision layer

| 层级 | 当前路径 | 作用 |
| --- | --- | --- |
| current freeze | `phase_reference_current_state_freeze_v1/` | 汇总当前 C1/C2/endpoint evidence 的状态冻结 |
| limited evidence | `phase_reference_limited_evidence_audit_v1/` | raw-curve coverage 和 shape audit |
| manual overlay | `phase_reference_manual_overlay_promotion_audit_v1/` | manual CEP/curve overlay 的晋升阻塞审计 |
| raw archive pointer | `raw_curve_archive_v1/` | 完整 raw curve 外部归档的 provenance 入口 |

`phase_reference_*` 是汇总/决策层。它们可以引用下游 case，但不能取代下游 manifest、hash、失败点或 unresolved 语义。冻结包记录生成时的路径和 hash；本次 namespace 整理不重写这些历史快照。

## Reading Order

1. 先读对应 case 的 `README.md` 和 `manifest.json`。
2. 再读 `AUDIT.md`、`decision.json` 和 `tables/claim_ledger.csv`。
3. 最后查看 figures 和派生表；不要把图面闭合、CEP midpoint 或 `s_shape_present` 自动解释为 physical promotion。
