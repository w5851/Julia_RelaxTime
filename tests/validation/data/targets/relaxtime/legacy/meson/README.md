# Legacy meson validation targets

This directory stores Fortran exported meson-level reference data
for Mott and meson-mass validation.

Expected CSV files (provided by legacy export pipeline):

- `legacy_meson_scan_fortran_muB0_v1.csv`
- `legacy_meson_scan_fortran_muB600_v1.csv`

See export spec:

- `docs/dev/active/2026-03-20_legacy_meson_export_spec.md`

Quality governance note:

- Rows with `solver_status=excluded_low_quality_nonconverged` are retained for audit,
  but should be excluded from strict numeric equality validation because the original
  legacy root diagnostics indicate non-converged / low-quality solve state.
