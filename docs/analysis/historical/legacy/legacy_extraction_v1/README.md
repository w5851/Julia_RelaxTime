# Legacy extraction and cross-check v1

## Scope

This artifact records a lightweight extraction from the eight directories that
were audited under `D:\Desktop\legacy` on 2026-08-23. It is an audit record,
not a production validation dataset and not a source-code import.

The purpose is to retain the useful conclusions after the external legacy
trees are removed:

- identify which legacy material is already represented by current targets,
  provenance records, formulas, or APIs;
- preserve the small set of numerical and implementation facts that affect
  current engineering decisions;
- record why the remaining source, binaries, and large scan outputs are not
  admitted to the main project.

## Audit provenance

- Extraction baseline commit: `f6146a8415c7e93ef0a023e1cf6550790791eded`
- Current main checked commit: `1ccf29310fb20c30bcd154f0b4966e25a7565225`
- Drift check: the current-main changes after the extraction baseline are
  phase-reference consumer adapters and the runtime switch only; they do not
  change the magnetic, Mott, RMFT, RelaxTime legacy target, or legacy inventory
  surfaces audited here. The observations below therefore remain current as of
  the checked commit.
- External source root: `D:\Desktop\legacy`
- Inventory: `3,184` files, `717,610,994` bytes across eight directories.
- Directory hashes and the hash construction rule are in
  `manifest.json` and `tables/legacy_extraction_index_v1.csv`.
- The tree hash is SHA-256 over UTF-8 lines sorted by relative path, where
  each line is `relative_path|file_sha256|file_size` and lines are joined by
  LF without an additional trailing line.

## What was extracted

The index and route audit retain only observations that can be traced to a
specific legacy artifact or to a current-project path:

- the archive snapshot is byte-identical to the gas-liquid and transport
  source trees at the file-content level;
- the Meson export contains 112 `exact_light_pole`, 40
  `exact_light_no_bracket`, and 456 `proxy_mass_from_threshold` rows; proxy
  rows are not formal meson-mass or Mott evidence;
- the Mott A/G/K replay has a maximum absolute difference of about `1.01e-7`
  and maximum relative difference about `1.64e-8`, while B0 has an imaginary
  difference of about `3.61e-10` and a real-part maximum relative difference
  about `0.902%`;
- the magnetic legacy route contains a charge-variable mismatch in
  thermodynamic postprocessing, NaNs in a fixed-density output, inconsistent
  executable/source variants, and an encoding-dependent equal-entropy run;
- the RMFT T-rho output contains 19,773 fixed-width overflow rows out of
  20,001 rows; RelaxTime legacy runs include non-finite output, zero-byte
  outputs, and a singular-matrix stop;
- no legacy code, executable, full scan, or unmanifested output is admitted
  as a current production input by this artifact.

The magnetic external-source comparison is recorded separately in
[`magnetic_external_crosscheck_v1.md`](magnetic_external_crosscheck_v1.md) and
[`tables/legacy_magnetic_external_crosscheck_v1.csv`](tables/legacy_magnetic_external_crosscheck_v1.csv).
The short FixedMu replay and same-route Julia/Fortran comparison are recorded in
[`magnetic_external_crosscheck_v2.md`](magnetic_external_crosscheck_v2.md) and
[`tables/magnetic_external_fixedmu_replay_v2.csv`](tables/magnetic_external_fixedmu_replay_v2.csv).
The aligned state deltas are in
[`tables/magnetic_julia_fortran_same_route_v2.csv`](tables/magnetic_julia_fortran_same_route_v2.csv).
The later paper-route audit and `pnjl_cep` test-set handoff are recorded in
[`magnetic_reference_route_audit_v3.md`](magnetic_reference_route_audit_v3.md).
The current Mott integral audit statistics are recorded in the main validation
provenance file
[`tests/validation/data/provenance/relaxtime/evidence/relaxtime_mott_integral_external_crosscheck_v1.csv`](../../../../../tests/validation/data/provenance/relaxtime/evidence/relaxtime_mott_integral_external_crosscheck_v1.csv).
Neither record promotes the external source to acceptance: the magnetic
parameter/unit/ensemble/cutoff/branch gates and the Mott convention gate remain
open. The v2 replay supplies only same-kernel diagnostic states; the paper-route
audit now marks smooth-Landau as non-acceptance until the regularization is
resolved. The next external-source target is `pnjl_cep`, pending a source
URL/commit and machine-readable output contract. The original point list is recorded in
[`magnetic_external_validation_point_plan_v1.csv`](tables/magnetic_external_validation_point_plan_v1.csv).

## Current-project relationship

The current project already carries the reusable lightweight pieces in the
following places:

- Meson/Mott targets and tests: `tests/validation/data/targets/relaxtime/legacy/meson/`
  and `tests/validation/relaxtime/`;
- transport target/provenance families:
  `tests/validation/data/targets/relaxtime/legacy/transport/` and
  `tests/validation/data/provenance/relaxtime/`;
- current PNJL thermal and solver implementations:
  `src/models/pnjl_physics/`, `docs/reference/formula/models/pnjl/`, and
  `docs/reference/formula/models/pnjl_magnetic/`;
- current gas-liquid formula documentation:
  `docs/reference/formula/models/gas_liquid/GasLiquid_RMF_CoreEquations.md`.

The mapping and admission decisions are explicit in
`tables/legacy_current_project_absorption_v1.csv`. Existing target files are
not overwritten or regenerated by this extraction.

## Cleanup gate

Deletion is intentionally a separate operation. Before removing the external
trees, review this artifact, verify the recorded inventory against the live
paths, and confirm that no untracked legacy file is still needed. The planned
cleanup target is exactly the eight child directories listed in the index;
the main project and unrelated desktop data are outside the cleanup scope.

After deletion, the artifact remains as the provenance record. It does not
claim that the legacy executables remain reproducible or that diagnostic-only
results are formal validation.
