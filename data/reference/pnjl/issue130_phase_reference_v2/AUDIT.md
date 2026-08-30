# Issue #130 phase-reference v2 solver-free materialization audit

The author has accepted the `accepted` layer as the downstream analysis default. This does not upgrade interpolated rows to strict certification and does not change solver runtime selection.

- strict/render/accepted source values remain unchanged except for the accepted row-level decision metadata;
- `solver_called=false`, `reference_write=false`, and `runtime_consumption=false`;
- strict runtime and the explicit PNJL legacy fallback/rollback remain available.
