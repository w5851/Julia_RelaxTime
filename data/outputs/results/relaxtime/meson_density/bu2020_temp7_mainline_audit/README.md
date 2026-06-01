# BU2020/temp7 mainline meson-density audit

date: 2026-06-01

## Scope

- temp7 evidence: `D:\Desktop\Julia_RelaxTime\temp7\particles2020_bu_kpi_independent_audit`
- output CSV: `data\outputs\results\relaxtime\meson_density\bu2020_temp7_mainline_audit\bu2020_temp7_meson_density_audit.csv`
- state point: `T=90.0 MeV`, `mu_q=350.0 MeV`, `xi=0.0`
- flavor profile: `bu2020_mu_s_0p2`

## Interpretation

- `real_axis_mode=:pv_b0_eta0` is emitted as a separate branch from finite-eta smoothing.
- `density_policy=:strict_normal_domain` marks Bose-domain unsafe rows instead of silently skipping `omega <= mu_M`.
- `density_policy=:excitation_only_E_gt_mu` is an explicit diagnostic continuation, not a literature fact.
- no-anomalous subtraction is not enabled in this mainline audit output.

## Status Counts

- `ok`: `12`
- `unsafe_bose_domain`: `8`
