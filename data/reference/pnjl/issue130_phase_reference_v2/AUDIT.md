# Issue #130 phase-reference v2 solver-free materialization audit

- source candidate: `data/reference/pnjl/issue130_phase_reference_v1`
- calculation SHA: `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- strict rows: boundary=7162, crossover=1343, cep=93, spinodals=6886
- render rows: boundary=12537, crossover=3135, cep=161, spinodals=11989
- accepted rows: boundary=12537, crossover=3135, cep=161, spinodals=11989
- render spinodal is copied from the immutable derived build input because v1 render had no spinodal coordinate table.
- accepted rows retain `strict_certified`, `unresolved` and `interpolated_noncertified`; no status is upgraded.
- `solver_called=false`, `oracle_labels_consumed=false`, `reference_write=false`, `runtime_consumption=false`.

## Boundary

This package is a candidate for author-reviewed downstream research. It is not a strict convergence certificate and does not authorize runtime switching or PNJL legacy physical deletion.
