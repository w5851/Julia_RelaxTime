# Validation Layout

`tests/validation/` stores acceptance-style numerical validation and baseline checks.

## Data Taxonomy

Validation data under `tests/validation/data/` is split into three semantic layers:

- `raw_long/`: long-form digitized tables, scan outputs, and large source datasets.
- `targets/`: lightweight CSV files consumed directly by `test_*.jl` acceptance checks.
- `provenance/`: disputed points, replacement evidence, and legacy consensus audit records.

Tests should default to `targets/`. Direct reads from `raw_long/` are reserved for explicit raw-compare workflows.

## Path Helpers

Shared path helpers live in `tests/validation/common/data_paths.jl`:

- `validation_raw_long_path(...)`
- `validation_targets_path(...)`
- `validation_provenance_path(...)`

New validation code should use these helpers rather than hardcoding `tests/validation/data/...` paths.

## Test Discovery

`tests/validation/runtests.jl` recursively includes nested `test_*.jl` files while preserving the top-level `PNJL Baselines` and `RelaxTime Baselines` testsets.

This allows future subdivision such as:

- `pnjl/literature/`
- `pnjl/reference/`
- `relaxtime/literature/`
- `relaxtime/legacy/`
- `relaxtime/advanced/`

without changing the entrypoint again.

## Family Semantics

- `tests/validation/pnjl/literature/`: PNJL acceptance checks anchored to literature text, figures, or thesis claims.
- `tests/validation/pnjl/reference/`: PNJL acceptance checks anchored to repository reference artifacts after extracting lightweight fixed-point targets.
- `tests/validation/relaxtime/legacy/transport/`: RelaxTime legacy transport guardrails anchored to legacy Fortran exports, with one quantity family per test file.

Matching target/provenance directories should preserve the same semantic split:

- `targets/pnjl/reference/`: lightweight crossover or phase fixed-point targets consumed directly by tests.
- `targets/relaxtime/legacy/transport/`: lightweight acceptance targets for legacy transport quantities such as `tau`, `eta_over_s`, `sigma_t`, `kappa_*`, `lambda`, and stable ratio diagnostics.
- `provenance/.../evidence/`: source file paths, column semantics, and comparison notes that justify how targets were extracted.

Practical rule:

- tests read `targets/`
- large scan outputs or reference curves stay out of default test input
- provenance records explain where targets came from and why a quantity is or is not in acceptance

## Admission Rules

New validation families must satisfy an explicit admission split before landing:

- `acceptance`: lightweight machine-readable targets that are stable enough to gate default validation.
- `evidence`: source audits, disputed points, rejected legacy anchors, branch-divergence notes, or empty-source inventories that justify why something is or is not in acceptance.
- `analysis`: large scans, decomposition tables, and exploratory outputs that are useful for diagnosis but should not become default test inputs.

Acceptance targets should satisfy all of the following:

- target rows are lightweight fixed points or similarly small representative subsets, not raw full scans
- field semantics map directly onto a stable Julia public or validation-facing quantity
- tolerance policy is explicit in the target file or family documentation
- source provenance is recorded under `provenance/`
- points are branch-unique or otherwise thermodynamically unambiguous for the intended comparison

Evidence-only records should be used when any of the following holds:

- legacy sources disagree and no source-consistency gate has been passed
- compared points are branch-ambiguous, metastable, or known to sit on different equilibrium branches
- checked-in legacy/reference exports are empty or schema-only
- the quantity currently exists only as a diagnostic or exploratory output rather than a stable acceptance contract

Analysis-only inputs should be used when any of the following holds:

- default validation would need the full raw scan or decomposition table to make sense of the result
- the comparison still depends on one-off scripts or manual interpretation
- the quantity is useful for research/debugging but not yet mature enough for fixed-point guardrails

Practical guardrail:

- a disputed or legacy point should not enter `targets/` just because it is numerically available
- admission requires source consistency, branch clarity, and a stable Julia-side quantity definition

## PNJL Crossover Mini Family

The PNJL crossover reference area currently has three distinct roles and they should not be merged conceptually:

- `pnjl_crossover_reference_targets_v1.csv`: repository reference acceptance targets extracted from the Julia-side reference artifact
- `pnjl_crossover_legacy_dual_source_targets_v1.csv`: legacy dual-source guardrail targets used only after shared fixed points have been judged source-consistent and branch-safe
- `pnjl_crossover_legacy_source_consistency_targets_v1.csv`: source-consistency gate that checks whether legacy Fortran and legacy C++ agree with each other before those points are treated as acceptance-worthy guardrails

For this mini family, the admission order is fixed:

1. legacy source points first land as provenance/evidence
2. shared points become source-consistency targets only after both sources are machine-readable and semantically aligned
3. only source-consistent and branch-unique fixed points may be promoted into the dual-source guardrail acceptance set

Current explicit exclusions for crossover acceptance:

- raw legacy scan curves are not default test inputs
- single-source legacy points are not promoted by default if a second source is missing
- high-temperature or low-density branch-ambiguous points are not valid acceptance anchors just because an old code exported them
- legacy points that conflict because of branch mismatch belong under `provenance/.../evidence/`, not under `targets/`

## Current Naming Conventions

- Raw long tables: `*_digitized_longtable_v1.csv`, `*_scan_longtable_v1.csv`, `*_legacy_longtable_v1.csv`
- Acceptance targets: `*_literature_targets_v1.csv`, `*_legacy_consensus_targets_v1.csv`, `*_legacy_state_targets_v1.csv`
- Provenance records: `*_disputed_points_v1.csv`, `*_replacement_points_v1.csv`, `*_evidence_index_v1.csv`

Field guidance for current families:

- RelaxTime legacy transport scalar targets use: `record_id,T_MeV,muB_MeV,field,expected_value,rtol,source`
- PNJL crossover reference targets use: `target_id,observable,xi,mu_MeV,expected_T_MeV,lower_T_MeV,upper_T_MeV,method,variable,source`