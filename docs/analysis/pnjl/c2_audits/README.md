# C2 Audit Evidence

This group contains the two cross-layer C2 audit packages:

- `c2_convergence_audit_v1/`: the C0/C1/C2 convergence and classification comparison audit.
- `c2_blocking_audit_v2/`: the current C1/C2 blocking, failure-reason, and unresolved-state audit.

They answer different audit questions and keep separate source runs, calculation SHAs, tables, figures, and verdicts. This directory groups them physically without merging their provenance or changing their diagnostic-only boundaries.

The C2 surface views remain under `../c2_surface_views/`; targeted review and feasibility plans remain separate follow-up groups.
