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
- Original upstream branch policy: `trho_reverse_rho=true`, `trho_seed_policy=temperature_grouped_rho_continuity`
- Targeted branch repair: three known branch-unstable points were recomputed on 2026-06-20 with `trho_branch_policy=pressure_max_all_attempts_multiseed`.

## Production Outputs

- CSV: `data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.csv`
- Run command: `data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/run.command.txt`
- Remote manifest: `data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/remote_run_manifest.json`
- Figure directory: `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/`
- SVG: `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.svg`
- PNG: `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.png`
- Plot manifest: `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/plot_manifest.json`

## Rendering

The original figure-side SVG was a rendering issue: it used `muq_MeV`
as the heatmap x-axis, while this production case is scanned on the regular
`rho_target` axis. After the 2026-06-20 branch repair described below, both
CSV and figures were updated.

Corrected plot-only commands:

```text
python scripts/analysis/relaxtime/render_combined_meson_density_fig3_like.py --csv data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.csv --out data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.svg --manifest data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/plot_manifest.json --field kpi_ratio --x-field rho_target --x-label "rho/rho0" --x-unit "" --kind trho_asymmetric_t_rho_heatmap_log_svg --color-scale log --dpi 260 --title "trho_asymmetric K+/pi+ T-rho heatmap (log color scale)"
python scripts/analysis/relaxtime/render_combined_meson_density_fig3_like.py --csv data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.csv --out data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.png --manifest data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/plot_manifest.json --field kpi_ratio --x-field rho_target --x-label "rho/rho0" --x-unit "" --kind trho_asymmetric_t_rho_heatmap_log_png --color-scale log --dpi 260 --title "trho_asymmetric K+/pi+ T-rho heatmap (log color scale)"
```

The corrected heatmap uses `rho_target` for x, `T_MeV` for y, and
`kpi_ratio` for color with a logarithmic color scale. The log scale is used
only for visualization and does not change the stored data.

## Targeted Branch Repair

On 2026-06-20, this artifact was updated after the upstream
`FixedAsymmetricRho` multi-root handling was fixed. The original production
scan selected a lower-pressure equilibrium branch at three isolated points.
Those points were recomputed locally with the same v1 numerical parameters and
the current upstream policy `pressure_max_all_attempts_multiseed`; only the
corresponding 12 CSV rows were replaced.

Repair evidence is stored under:

- `data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/branch_repair_20260620/`

Repaired points:

| point | old issue | repaired `mu_u,mu_d,mu_s` MeV | repaired `phase_shift_current K/pi` | repaired `GBU K/pi` |
| --- | --- | ---: | ---: | ---: |
| `T=120 MeV, rho=0.35` | lower-pressure branch, extreme ratio | `290.818, 299.261, 31.377` | `82.607025` | `129.863384` |
| `T=130 MeV, rho=0.80` | lower-pressure branch, extreme ratio | `288.492, 298.823, 16.971` | `88.127569` | `137.065471` |
| `T=210 MeV, rho=0.90` | weak branch/order-sensitive point in old artifact | `46.552, 53.019, 0.121` | `12.372666` | `1.144906` |

The rest of the grid was not recomputed in this repair pass. New branch
diagnostic fields are therefore populated for repaired rows and may be empty for
original rows.

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

Production rows: `880` total, with `ok=661` and `unsafe_bose_domain=219`.

By regime:

- `stable`: `ok=208`, `unsafe_bose_domain=12`
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
