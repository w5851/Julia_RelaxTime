# Combined meson-density T-mu scan

date: 2026-06-01

## Scope

- path strategy: `tmu`
- density regimes: `stable`, `strict_bw_stage1`, `phase_shift_current`, `phase_shift_gbu_reference`
- fixed-mu path values: `mu_q=0.000000,50.000000,100.000000,150.000000,200.000000,250.000000,300.000000 MeV`
- temperature range: `120.0:10.0:240.0 MeV`
- flavor profile: `default`
- meson profile: `default`

## Outputs

- CSV: `data\outputs\results\relaxtime\meson_density\fig3_like_tmu_equal_flavor_folded_highres\combined_meson_density_scan.csv`
- SVG: `data\outputs\results\relaxtime\meson_density\fig3_like_tmu_equal_flavor_folded_highres\combined_meson_density_scan.svg`

## Regime Definitions

- `stable`: stable particle limit from `Models.solve_meson_density_from_meson_point`.
- `strict_bw_stage1`: reduced strict-BW density from `Models.solve_strict_bw_meson_density_from_meson_point`.
- `phase_shift_current`: phase-shift BU density with weight `delta`.
- `phase_shift_gbu_reference`: generalized BU density with weight `delta - 0.5sin(2delta)`.

## Policies

- `real_axis_mode=pv_b0_eta0`
- `phase_convention=arg_inverse_propagator`
- `phase_display=fold_0_pi`
- `density_policy=strict_normal_domain`
- `noanom_policy=none`
- `gamma_zero_tol=1.0e-6`
- This entrypoint is a Bridge-style composition of scan path and density-regime strategy.

## Status Counts

- `ok`: `364`
