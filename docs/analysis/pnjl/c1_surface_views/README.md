# C1 Surface Views

This group contains two diagnostic projections of the same C1 source run:

- `pnjl_c1_mu_xi_T_phase_surfaces_diagnostic_v2/`: `(mu, xi, T)` view with the CEP temperature-bracket figure.
- `pnjl_c1_xi_t_mu_phase_surfaces_diagnostic_v1/`: `(xi, T, mu)` view.

Both views use source workflow run `31762201725` and the same C1 input artifact. They remain separate because their axis order, plot manifest schema, and figure semantics differ. The packages are `diagnostic_only`; unresolved grid records and ambiguous CEP rows are not promoted to physical claims.

The move only changes repository paths. Input tables, derived tables, manifests, and PNG bytes were not recomputed or redrawn.
