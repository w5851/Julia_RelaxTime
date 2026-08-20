# C2 Follow-up Evidence

This directory groups the follow-up packages that consume or extend the C2 audit evidence. They answer different questions and remain separate provenance nodes:

- `c2_targeted_manual_review_v1/`: solver-free review of nine targeted rho-mu curves plus the three CEP midpoint plans; author review is still required.
- `c2_limited_feasibility_v1/`: the frozen 17-bracket CEP input contract for the limited-feasibility workflow; it is not a numerical result.
- `c2_cep_xi05_high_side_extension_v1/`: the diagnostic `xi=0.5` high-side temperature extension plan.
- `c2_manual_bisection_v1/`: the fixed low/mid/high temperature plan for the three unresolved CEP brackets.

Each package keeps its own README, manifest, tables, figures, hashes, and diagnostic boundary. The physical grouping does not merge evidence, rewrite source provenance, call a solver, or promote a production/reference result. `phase_reference_*` packages remain at the PNJL root as decision and freeze layers.
