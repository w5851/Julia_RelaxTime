# trho_asymmetric K+/pi+ Meson Density Scan v1

Verdict: `production-grade`

This directory contains the formal `FixedAsymmetricRho` density-constrained
`K+ / pi+` meson-density production artifact.

## Scope

- Path strategy: `trho_asymmetric`
- Equilibrium source: `FixedAsymmetricRho`
- Meson chemical profile: `asymmetric_kplus_over_piplus_signed`
- Constraint targets: `rho_u / rho_d = 0.876`, `rho_s = 0 fm^-3`
- Grid: `T = 120:10:220 MeV`, `rho/rho0 = 0.05:0.05:1.00`
- Regimes: `stable`, `strict_bw_stage1`, `phase_shift_current`, `phase_shift_gbu_reference`
- Phase-shift policy: `density_policy=x_min_cut`, `bose_x_min=1e-2`
- Policy scope: `x_min_cut` applies only to BU/GBU phase-shift regimes; stable/BW keep strict Bose-domain diagnostics.
- Upstream branch policy: `trho_reverse_rho=true`, `trho_seed_policy=temperature_grouped_rho_continuity`

## Production Outputs

- CSV: `data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.csv`
- Run command: `data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/run.command.txt`
- Remote manifest: `data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/remote_run_manifest.json`
- Figure directory: `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/`
- SVG: `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.svg`
- PNG: `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.png`
- Plot manifest: `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/plot_manifest.json`

## Rendering

The original figure-side SVG was a rendering-only issue: it used `muq_MeV`
as the heatmap x-axis, while this production case is scanned on the regular
`rho_target` axis. The CSV data were not recomputed.

Corrected plot-only commands:

```text
python scripts/analysis/relaxtime/render_combined_meson_density_fig3_like.py --csv data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.csv --out data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.svg --manifest data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/plot_manifest.json --field kpi_ratio --x-field rho_target --x-label "rho/rho0" --x-unit "" --kind trho_asymmetric_t_rho_heatmap_svg --dpi 260 --title "trho_asymmetric K+/pi+ T-rho heatmap"
python scripts/analysis/relaxtime/render_combined_meson_density_fig3_like.py --csv data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.csv --out data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.png --manifest data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/plot_manifest.json --field kpi_ratio --x-field rho_target --x-label "rho/rho0" --x-unit "" --kind trho_asymmetric_t_rho_heatmap_png --dpi 260 --title "trho_asymmetric K+/pi+ T-rho heatmap"
```

The corrected heatmap uses `rho_target` for x, `T_MeV` for y, and
`kpi_ratio` for color. A few bright cells are real large-ratio rows in the CSV,
not residual coordinate artifacts.

## Convergence Gate

The selected production parameters are `q_nodes=192`, `omega_nodes=192`,
`stable_q_nodes=768`, `p_num=8`, `t_num=4`, `max_iter=20`.

The final adjacent high-precision check is `custom128 -> custom192`:

| field | comparable ok rows | max relative diff | mean relative diff | max point |
| --- | ---: | ---: | ---: | --- |
| `n_pi` | 660 | 0.00930591 | 0.000841249 | `{'T_MeV': 120.0, 'rho_target': 0.4, 'regime': 'phase_shift_gbu_reference'}` |
| `n_K` | 660 | 0.00924145 | 9.65222e-05 | `{'T_MeV': 220.0, 'rho_target': 0.2, 'regime': 'strict_bw_stage1'}` |
| `kpi_ratio` | 660 | 0.00926615 | 0.000896029 | `{'T_MeV': 120.0, 'rho_target': 0.4, 'regime': 'phase_shift_gbu_reference'}` |

Upstream branch stability in the same check: `muq_MeV`, `m_pi_MeV`, and
`m_K_MeV` all have max absolute difference `0.0` on stable rows.

## Status Counts

Production rows: `880` total, with `ok=660` and `unsafe_bose_domain=220`.

By regime:

- `stable`: `ok=207`, `unsafe_bose_domain=13`
- `strict_bw_stage1`: `ok=13`, `unsafe_bose_domain=207`
- `phase_shift_current`: `ok=220`
- `phase_shift_gbu_reference`: `ok=220`

`unsafe_bose_domain` rows are retained as explicit strict-domain diagnostics and
must not be interpreted as finite stable/BW densities.

## Files

- `combined_meson_density_scan.csv`
- `README.md`
- `PRODUCTION_AUDIT.md`
- `manifest.json`
- `remote_run_manifest.json`
- `run.command.txt`
- `run.exitcode`
- `run.stdout.log`
- `run.stderr.log`
- `convergence/convergence_summary.json`
