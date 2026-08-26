# first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2 (mode_b_fixed_T_sparse_muB)

This is the versioned Issue #130 RS candidate result imported from the
author-reviewed numerical aggregate. It is not the default runtime result.

## Status

- artifact status: `imported_candidate`
- numerical status: `diagnostic_only`
- runtime default switch: `false`
- legacy fallback: retained; explicit rollback is `--phase-reference-mode legacy`
- source calculation SHA: `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- source workflow head: `22874505877491754eed27519ad8a7b871c82571`

## Scope

- mode: `mode_b_fixed_T_sparse_muB` (fixed temperature with sparse mu_B values)
- aggregate scan rows: `909`
- selected source shards: `15`
- direct-coexistence policy: `xi=+/-0.003` retained for the mode-A `muB=900, alpha_T=1` anchor

## Provenance

- immutable aggregate: `aggregate_replay_20260826_v4`
- aggregate manifest SHA256: `42e88e028832311a95e45816eaf8ab1bb2dbbe8bb094ed6b59681c51d1cad754`
- post-repair audit: `post_repair_audit_20260826_v1`
- post-repair verdict: `post_repair_audit_pass_diagnostic_only`
- source numerical solver was used; this import and aggregate replay were solver-free

## Quality and limitations

Common `tau_u_ubar_ratio_high` warnings are retained in the scan and manifest.
They are diagnostic warnings, not silently filtered failures. The known
historical sidecar hash defect remains recorded; PR #269 fixed the producer
contract but did not rewrite immutable historical artifacts.

This case does not switch runtime consumers, delete `prod_v1`, or claim that
candidate/legacy numerical differences are a convergence result.
