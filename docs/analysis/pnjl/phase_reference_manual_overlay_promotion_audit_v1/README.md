# Phase-reference manual overlay promotion audit v1

This solver-free package evaluates whether C2 automatic evidence plus the nine-point and CEP author overlays can support a limited phase-reference promotion. It does not write `data/reference/**` and does not change C2 artifacts.

The current verdict is `blocked_manual_overlay_inconclusive`. The C2 grid contains 5424 unresolved rows; 4661 are geometry/interpolation-only by reason, while 763 remain unclassified, Stage-C, or classification-transition blockers. The nine oracle curves and three CEP brackets are therefore retained as explicit diagnostic overlays rather than silently replacing hybrid states.

See `decision.json`, `execution_log.md`, `tables/auto_gate_summary.csv`, `tables/grid_unresolved_by_reason.csv`, `tables/manual_nine_curve_decisions.csv`, `tables/manual_cep_decisions.csv`, and `tables/claim_ledger.csv`.
