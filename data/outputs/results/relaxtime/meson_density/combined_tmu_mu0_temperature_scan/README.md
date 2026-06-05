# Combined meson-density T-mu scan

date: 2026-06-05

## Scope

- path strategy: `tmu`
- density regimes: `stable`, `strict_bw_stage1`, `phase_shift_current`, `phase_shift_gbu_reference`
- fixed chemical potential: `mu_q=0.0 MeV`, `mu_B=0.0 MeV`
- temperature range: `120.0:10.0:240.0 MeV`
- flavor profile: `default`
- meson profile: `default`

## Outputs

- result directory: `data/outputs/results/relaxtime/meson_density/combined_tmu_mu0_temperature_scan`
- CSV: `data/outputs/results/relaxtime/meson_density/combined_tmu_mu0_temperature_scan/combined_meson_density_scan.csv`
- figure directory: `data/outputs/figures/relaxtime/meson_density/combined_tmu_mu0_temperature_scan`
- SVG: `data/outputs/figures/relaxtime/meson_density/combined_tmu_mu0_temperature_scan/combined_meson_density_scan.svg`
- PNG: `data/outputs/figures/relaxtime/meson_density/combined_tmu_mu0_temperature_scan/combined_meson_density_scan.png`
- plot manifest: `data/outputs/figures/relaxtime/meson_density/combined_tmu_mu0_temperature_scan/plot_manifest.json`

## Rendering


- PNG command: `python scripts/analysis/relaxtime/render_combined_meson_density_temperature_scan.py --csv data/outputs/results/relaxtime/meson_density/combined_tmu_mu0_temperature_scan/combined_meson_density_scan.csv --out data/outputs/figures/relaxtime/meson_density/combined_tmu_mu0_temperature_scan/combined_meson_density_scan.png --manifest data/outputs/figures/relaxtime/meson_density/combined_tmu_mu0_temperature_scan/plot_manifest.json --readme data/outputs/results/relaxtime/meson_density/combined_tmu_mu0_temperature_scan/README.md --dpi 260 --title "Combined meson-density mu_q=0 temperature scan"`
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
- `gamma_zero_tol=1.0e-6`
- This entrypoint is a Bridge-style composition of scan path and density-regime strategy.

## Status Counts

- `ok`: `52`
