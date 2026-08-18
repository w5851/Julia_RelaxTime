# Phase-reference manual overlay promotion audit v1

- Generated (UTC): `2026-08-17T13:25:23.209198+00:00`
- Command: `scripts/analysis/pnjl/build_phase_reference_manual_overlay_promotion_audit.py --c1-root D:/Desktop/Julia_RelaxTime_issue130_artifacts/c1_complete_acceptance_31762201725_attempt2 --c2-root D:/Desktop/Julia_RelaxTime_issue130_artifacts/stagec_density_v2_c2_20260813_run31862752226 --c2-audit-root D:/Desktop/Julia_RelaxTime/docs/analysis/pnjl/c2_blocking_audit_v2 --nine-review-root D:/Desktop/Julia_RelaxTime/docs/analysis/pnjl/c2_targeted_manual_review_v1 --manual-overlay D:/Desktop/Julia_RelaxTime_issue130_artifacts/cep_manual_bisection_audit_v2_31999149922/author_decision_overlay_v2.json --output-root D:/Desktop/Julia_RelaxTime/docs/analysis/pnjl/phase_reference_manual_overlay_promotion_audit_v1`
- Calculation SHA: `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- Solver called: `false`
- Reference write: `false`
- Input policy: C2 automatic artifacts plus explicit author-review overlays; no C2 row was overwritten.
- Result: `blocked_manual_overlay_inconclusive`
- Reason: unresolved grid rows include non-geometry categories; nine oracle acceptances and CEP overlay do not rewrite C2 hybrid labels.
