# Freeze-out Meson Mass Case: default_baseline_freezeout_xi0_loggrid_1to200_n30

Current role:

- formal freezeout meson-mass path output
- xi = 0.0
- sqrt(s_NN) sampled on a log10-uniform grid over [1, 200] GeV with 30 points

Known path-validity note:

- The formal production CSV is intentionally kept unchanged on `[1, 200] GeV`.
- For the `default` freezeout baseline profile, the first two sampled points produce nonpositive `T_MeV` and are therefore retained as explicit failure rows in the CSV.
- The default-profile `T=0` threshold is approximately `sqrt(s_NN) = 1.411898 GeV`.
- For future physical freezeout scans, use `sqrt(s_NN) >= 1.45 GeV` for this profile family.

Plotting policy:

- Figures in `figures/` exclude rows with `equilibrium_converged != true` or `T_MeV <= 0`.
- Excluded points in this case: `2`.

Files:

- source CSV: `data/outputs/results/relaxtime/meson_mass/path_scan/freezeout/default_baseline_freezeout_xi0_loggrid_1to200_n30.csv`
- valid-points CSV: `data/outputs/results/relaxtime/meson_mass/path_scan/freezeout/default_baseline_freezeout_xi0_loggrid_1to200_n30.valid_points.csv`
- provenance: `data/outputs/results/relaxtime/meson_mass/path_scan/freezeout/default_baseline_freezeout_xi0_loggrid_1to200_n30.provenance.json`
- path overview: `data/outputs/results/relaxtime/meson_mass/path_scan/freezeout/figures/default_baseline_freezeout_xi0_loggrid_1to200_n30_path_overview.png`
- pseudoscalar masses: `data/outputs/results/relaxtime/meson_mass/path_scan/freezeout/figures/default_baseline_freezeout_xi0_loggrid_1to200_n30_pseudoscalar_masses.png`
- scalar masses: `data/outputs/results/relaxtime/meson_mass/path_scan/freezeout/figures/default_baseline_freezeout_xi0_loggrid_1to200_n30_scalar_masses.png`

Quick summary:

- total points: `30`
- plotted valid points: `28`
- excluded invalid points: `2`
- first plotted point: `sqrt(s_NN) = 1.441092 GeV`
- last plotted point: `sqrt(s_NN) = 200.000000 GeV`
