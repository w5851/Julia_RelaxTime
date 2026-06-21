# Combined meson-density T-mu scan

date: 2026-06-20

## Scope

- path strategy: `trho_asymmetric`
- density regimes: `stable`, `strict_bw_stage1`, `phase_shift_current`, `phase_shift_gbu_reference`
- FixedAsymmetricRho rho targets: `0.900000`
- FixedAsymmetricRho branch policy: temperature-grouped scan with continuation seed included in an all-attempt MultiSeed candidate pool; selected by pressure max under constraints.
- FixedAsymmetricRho scan order: `trho_reverse_rho=true`.
- asymmetry targets: `rho_u/rho_d=0.876`, `rho_s=0.0 fm^-3`
- production status: `trho_asymmetric` may be used for formal artifacts only after an explicit convergence gate and production audit; ad hoc runs remain diagnostic evidence.
- temperature range: `210.0:10.0:210.0 MeV`
- flavor profile: `default`
- meson profile: `asymmetric_kplus_over_piplus_signed`

## Outputs

- result directory: `data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/branch_repair_20260620/T210_rho090`
- CSV: `data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/branch_repair_20260620/T210_rho090/combined_meson_density_scan.csv`

## Regime Definitions

- `stable`: stable particle limit from `Models.solve_meson_density_from_meson_point`.
- `strict_bw_stage1`: reduced strict-BW density from `Models.solve_strict_bw_meson_density_from_meson_point`.
- `phase_shift_current`: phase-shift BU density with weight `delta`.
- `phase_shift_gbu_reference`: generalized BU density with weight `delta - 0.5sin(2delta)`.

## Policies

- `real_axis_mode=pv_b0_eta0`
- `phase_convention=arg_inverse_propagator`
- `phase_display=fold_0_pi`
- `density_policy=x_min_cut`
- `bose_x_min=0.01`
- `density_policy` and `bose_x_min` apply only to `phase_shift_current` and `phase_shift_gbu_reference`; `stable` and `strict_bw_stage*` keep strict Bose-domain diagnostics.
- `noanom_policy=none`
- `strict_bw_omega_min=0.05`
- `gamma_zero_tol=1.0e-12`
- This entrypoint is a Bridge-style composition of scan path and density-regime strategy.

## Status Counts

- `ok`: `3`
- `unsafe_bose_domain`: `1`
