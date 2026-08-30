# RS Candidate Result Import Audit (mode_b_fixed_T_sparse_muB)

## Verdict

`imported_candidate_diagnostic_only`

## Input gates

- calculation SHA: `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- workflow head: `22874505877491754eed27519ad8a7b871c82571`
- aggregate manifest: `42e88e028832311a95e45816eaf8ab1bb2dbbe8bb094ed6b59681c51d1cad754`
- post-repair audit manifest: `587f79db64476f777b25ca5ee234a43b5a31e17bea375864d9e13f0226d1d6f0`
- post-repair audit verdict: `post_repair_audit_pass_diagnostic_only`
- scan rows: `909`
- diagnostic rows: `38178`
- scan duplicate keys: `0`
- diagnostic duplicate keys: `0`
- solver/failed rows: `0` in selected source shards

## Provenance semantics

- `source_solver_called=true`: the selected numerical Actions produced the source data.
- `aggregate_replay_solver_called=false`: aggregate replay did not call the solver.
- `import_solver_called=false`: this repository import only copied and derived sidecars.
- `production_write=false`: no runtime/default production switch occurred.

## Preserved boundaries

- old `prod_v1` result trees are not modified;
- legacy phase-reference fallback and explicit rollback remain available;
- candidate/legacy differences remain input-difference diagnostics, not numerical parity;
- quality warnings and historical sidecar hash defects are preserved in `convergence/`.
