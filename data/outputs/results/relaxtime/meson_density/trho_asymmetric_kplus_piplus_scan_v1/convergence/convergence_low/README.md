# Combined meson-density T-mu scan

date: 2026-06-06

## Scope

- path strategy: `trho_asymmetric`
- density regimes: `stable`, `strict_bw_stage1`, `phase_shift_current`, `phase_shift_gbu_reference`
- FixedAsymmetricRho rho targets: `0.050000,0.100000,0.150000,0.200000,0.250000,0.300000,0.350000,0.400000,0.450000,0.500000,0.550000,0.600000,0.650000,0.700000,0.750000,0.800000,0.850000,0.900000,0.950000,1.000000`
- asymmetry targets: `rho_u/rho_d=0.876`, `rho_s=0.0 fm^-3`
- smoke-only status: this path is intended for diagnostic integration, not formal high-precision production.
- temperature range: `120.0:10.0:220.0 MeV`
- flavor profile: `default`
- meson profile: `asymmetric_kplus_over_piplus_signed`

## Outputs

- result directory: `data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/convergence/convergence_low`
- CSV: `data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/convergence/convergence_low/combined_meson_density_scan.csv`
- figure directory: `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/convergence/convergence_low`
- SVG: `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/convergence/convergence_low/combined_meson_density_scan.svg`
- plot manifest: `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/convergence/convergence_low/plot_manifest.json`

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
- `strict_bw_omega_min=0.05`
- `gamma_zero_tol=1.0e-12`
- This entrypoint is a Bridge-style composition of scan path and density-regime strategy.

## Status Counts

- `failed`: `221`
- `ok`: `245`
- `unsafe_bose_domain`: `414`
