# Production Audit: trho_asymmetric K+/pi+ Meson Density Scan v1

Verdict: `production-grade`

Date: 2026-06-08

## Production Case

- Case slug: `trho_asymmetric_kplus_piplus_scan_v1`
- Result path: `data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/`
- Figure path: `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/`
- Remote workflow: `.github/workflows/relaxtime-meson-density-production.yml`
- Production run: `27115936651`
- Source commit: `d2c5cba82702c18e55b2d398873145141a15a104`

## Physics Scope

- Equilibrium source: `FixedAsymmetricRho`
- Constraint targets: `sum(rho_i)/(3 rho0)=rho_target`, `rho_u/rho_d=0.876`, `rho_s=0 fm^-3`
- Charged profile: `asymmetric_kplus_over_piplus_signed`
- Meson channels: `pi+`, `K+`
- Grid: `T=120:10:220 MeV`, `rho/rho0=0.05:0.05:1.00`
- Density regimes: `stable`, `strict_bw_stage1`, `phase_shift_current`, `phase_shift_gbu_reference`

## Policies

- Upstream branch policy: `trho_reverse_rho=true`, `trho_seed_policy=temperature_grouped_rho_continuity`
- Phase policy: `real_axis_mode=pv_b0_eta0`, `phase_convention=arg_inverse_propagator`, `phase_display=fold_0_pi`
- BU/GBU density policy: `density_policy=x_min_cut`, `bose_x_min=1e-2`
- Scope of x-min policy: `phase_shift_current` and `phase_shift_gbu_reference` only
- Stable/BW policy: strict Bose-domain diagnostics; unsafe rows are emitted as `unsafe_bose_domain` / `NaN` and do not block BU/GBU production.

## Command Log

Production command is recorded in `run.command.txt`.

Selected production parameters:

- `p_num=8`, `t_num=4`, `max_iter=20`
- `stable_q_nodes=768`
- `qmax=4.0 fm^-1`, `q_nodes=192`
- `omega_min=0.05 fm^-1`, `omega_max=4.0 fm^-1`, `omega_nodes=192`
- `gamma_zero_tol=1e-12`

## Convergence Matrix

Remote runs used for convergence and production:

| tag | run id | role |
| --- | --- | --- |
| low | `27114834622` | coarse convergence |
| mid | `27114834644` | mid convergence |
| high36 | `27114834615` | initial high convergence |
| custom64 | `27115127952` | diagnostic refinement |
| custom96 | `27115242352` | diagnostic refinement |
| custom128 | `27115372364` | adjacent high reference |
| custom192 | `27115609958` | selected production precision reference |
| production192 | `27115936651` | formal production run |

Final adjacent check: `custom128 -> custom192`.

| field | max relative diff | mean relative diff | max point |
| --- | ---: | ---: | --- |
| `n_pi` | 0.00930591 | 0.000841249 | `{'T_MeV': 120.0, 'rho_target': 0.4, 'regime': 'phase_shift_gbu_reference'}` |
| `n_K` | 0.00924145 | 9.65222e-05 | `{'T_MeV': 220.0, 'rho_target': 0.2, 'regime': 'strict_bw_stage1'}` |
| `kpi_ratio` | 0.00926615 | 0.000896029 | `{'T_MeV': 120.0, 'rho_target': 0.4, 'regime': 'phase_shift_gbu_reference'}` |

Upstream branch stability in the final adjacent check:

- `muq_MeV_max_abs = 0.0`
- `m_pi_MeV_max_abs = 0.0`
- `m_K_MeV_max_abs = 0.0`

## Data Outputs

- `combined_meson_density_scan.csv`
- `remote_run_manifest.json`
- `run.command.txt`
- `run.exitcode`
- `run.stdout.log`
- `run.stderr.log`
- `convergence/convergence_summary.json`

## Figure Outputs

- `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.svg`
- `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.png`
- `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/plot_manifest.json`

Figure-side correction note:

- The first SVG quicklook used `muq_MeV` as the heatmap x-axis. For
  `trho_asymmetric`, `muq_MeV` is an equilibrium diagnostic field and is not
  the regular scan axis.
- The corrected SVG and PNG were rendered locally from the existing production
  CSV using `--x-field rho_target --x-label "rho/rho0" --x-unit ""` and
  `--color-scale log`.
- No density data, convergence evidence, or production CSV rows were recomputed
  for this plot-only correction.
- The largest bright cells are present in the source CSV: e.g.
  `phase_shift_gbu_reference` has `kpi_ratio=11481.086618` at
  `T=130 MeV, rho_target=0.8`, and `kpi_ratio=10833.443521` at
  `T=120 MeV, rho_target=0.35`.
- Those extreme ratios come from the quotient itself: at the two GBU points
  `n_pi=0.004475, n_K=51.382839` and `n_pi=0.004690, n_K=50.805504`,
  respectively. The same coordinates also carry nonzero `unsafe_bose_count`
  under the explicit BU/GBU `x_min_cut` policy, while strict stable/BW either
  remain unsafe or are only comparison diagnostics.
- The logarithmic color scale is a figure-only display policy chosen to keep
  the full finite range visible without clipping the extreme cells.

## Validation Commands And Results

Before the production run:

- `julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_combined_meson_density_scan_smoke.jl"; include("tests/integration/runtests.jl")'` passed.
- `julia --project=. scripts/dev/check_script_entrypoints.jl` passed.
- `julia --project=. scripts/dev/check_relaxtime_script_governance.jl` passed.
- `julia --project=. scripts/dev/check_docs_consistency.jl` passed.
- `julia --project=. scripts/dev/check_data_output_path_guard.jl` passed.
- `git diff --check` passed.

Slow local CLI smoke was attempted but timed out on this machine; remote GitHub
Actions runs listed above are the authoritative long-run evidence.

## Known Limitations And Residual Risks

- Stable and strict BW contain strict-domain `unsafe_bose_domain` rows under the selected `K+ / pi+` profile; those rows are not finite density values.
- BW remains a comparison regime and should not be used alone as the formal physical conclusion.
- The `x_min_cut` treatment is an explicit project policy for BU/GBU phase-shift regimes, not a unique prescription from the literature.
- The formal conclusion is scoped to this grid, profile, and policy; changing charged channel, asymmetry targets, or `bose_x_min` requires a new convergence gate.
