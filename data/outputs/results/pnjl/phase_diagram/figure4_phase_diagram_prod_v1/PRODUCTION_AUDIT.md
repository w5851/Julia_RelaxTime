# Production Audit: figure4_phase_diagram_prod_v1

Date: 2026-07-05

## Production Case

- case slug: `figure4_phase_diagram_prod_v1`
- result root: `data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/`
- figure root: `data/outputs/figures/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/`
- formal reference tag: `figure4_phase_diagram_prod_v1_c1_p24t8`
- formal reference verdict: `production-grade`

## Physics Scope

- PNJL phase reference for anisotropy scan `xi=-0.5:0.05:0.5`.
- Included observables: first-order boundary, CEP, crossover, spinodal branches.
- Default formal figure: two-panel `T-mu_B` and `T-rho` plot.
- Default plotted subset: `xi=-0.5,-0.25,0.0,0.25,0.5`.
- Default spinodal plotting policy: data included in production assets; visual layer disabled by default.

## Non-Goals

- No solver-physics rewrite.
- No transport-curve fitting or phase-boundary inference from transport observables.
- No paper-project edits, captions, claim tables, manuscript text, or LaTeX updates.
- No claim that GitHub Actions success alone is sufficient; production status is based on convergence and audit evidence.

## Command Log

GitHub Actions production/convergence runs:

- C1 full-grid candidate: https://github.com/w5851/Julia_RelaxTime/actions/runs/28736257287
- C2 refined anchors: https://github.com/w5851/Julia_RelaxTime/actions/runs/28741709519
- head SHA: `aca739f6bac4565dddef888abf9588b8c0fa583f`

Convergence comparison:

```powershell
python scripts/analysis/pnjl/compare_phase_reference_convergence.py --candidate-root <C1 artifact root> --candidate-tag figure4_phase_diagram_prod_v1_c1_p24t8 --reference-root <C2 artifact root> --reference-tag figure4_phase_diagram_prod_v1_c2_p32t12_anchors --xi-values -0.5,-0.25,0.0,0.25,0.5 --out-dir data/outputs/results/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/convergence/c1_vs_c2_anchor_comparison
```

Formal plot generation:

```powershell
python scripts/analysis/pnjl/build_figure4_phase_diagram_assets.py
```

Spinodal option smoke test:

```powershell
python scripts/analysis/pnjl/build_figure4_phase_diagram_assets.py --include-spinodal --result-dir <temp>/results --figure-dir <temp>/figures --formats png --dpi 80
```

## Convergence Matrix

| tier | coverage | numeric policy | rows | result |
|---|---|---|---:|---|
| C1 full-grid | `xi=-0.5:0.05:0.5` | `p_num=24`, `t_num=8`, `iterations=80` | boundary 304, CEP 21, crossover 336, spinodal 304 | accepted source |
| C2 refined anchors | `xi=-0.5,-0.25,0.0,0.25,0.5` | `p_num=32`, `t_num=12`, `iterations=100` | boundary 72, CEP 5, crossover 80, spinodal 72 | convergence check |

Key C1-vs-C2 maximum absolute differences:

| artifact | metric | max abs diff |
|---|---|---:|
| boundary | `mu_transition_MeV` | 0.2613591590 |
| CEP | `T_CEP_MeV` | 0.0292968750 |
| CEP | `muB_CEP_MeV` | 0.0759663347 |
| crossover | `T_crossover_MeV` | 0.0304144661 |
| spinodals | `mu_spinodal_hadron_MeV` | 0.0716456155 |
| spinodals | `mu_spinodal_quark_MeV` | 1.4341450249 |

Full convergence evidence is in `convergence/P4-C_CONVERGENCE_SUMMARY.md`.

## Selected Production Parameters

- `p_num=24`
- `t_num=8`
- `iterations=80`
- `T_min=60`, `T_max=240`, `T_step=5`
- `rho_min=0.0`, `rho_max=4.0`, `rho_step=0.05`
- `crossover_n_mu=16`
- `crossover_mu_max=450`

## Data Outputs

- `reference/boundary_figure4_phase_diagram_prod_v1_c1_p24t8.csv`
- `reference/cep_figure4_phase_diagram_prod_v1_c1_p24t8.csv`
- `reference/crossover_figure4_phase_diagram_prod_v1_c1_p24t8.csv`
- `reference/crossover_figure4_phase_diagram_prod_v1_c1_p24t8.meta.json`
- `reference/phase_reference_figure4_phase_diagram_prod_v1_c1_p24t8_manifest.json`
- `reference/spinodals_figure4_phase_diagram_prod_v1_c1_p24t8.csv`
- `phase_reference_source_manifest.json`
- `figure_assets/figure4_phase_lines_TmuB.csv`
- `figure_assets/figure4_phase_lines_Trho.csv`
- `figure_assets/figure4_phase_plot_inputs_summary.json`

The `T-rho` figure asset includes plotting-only `CEP_visual_connector` rows for each `xi`: two coexistence connector rows and one crossover connector row. These rows make the plotted coexistence branches and crossover curve visually pass through CEP; they do not modify the formal reference CSV files.

## Figure Outputs

- `data/outputs/figures/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/figure4_phase_diagram_TmuB_Trho.png`
- `data/outputs/figures/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/figure4_phase_diagram_TmuB_Trho.pdf`
- `data/outputs/figures/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/plot_manifest.json`

## Validation Results

Passed:

- `python -m py_compile scripts/analysis/pnjl/build_figure4_phase_diagram_assets.py scripts/analysis/pnjl/compare_phase_reference_convergence.py scripts/pnjl/validate_dense_reference_artifact.py`
- `python scripts/analysis/pnjl/build_figure4_phase_diagram_assets.py`
- temporary `--include-spinodal` smoke generation
- `T-rho` CEP visual connector check: 63 connector rows generated, with 15 plotted in the default five-`xi` subset
- plot manifest hash check: 10 files verified
- PNG dimension check: `2244x1156`, `RGBA`
- visual inspection: default figure is readable; axes, legend, panel labels, default no-spinodal density, and `T-rho` CEP line connectivity are acceptable for the current asset review stage
- `julia --project=. scripts/dev/check_docs_consistency.jl`
- `julia --project=. scripts/dev/check_script_entrypoints.jl`
- `julia --project=. scripts/dev/check_data_output_path_guard.jl`
- `git diff --check`

Known unrelated governance failure:

- `julia --project=. scripts/dev/check_active_docs_governance.jl` fails because two pre-existing active docs are stale: `2026-04-30_介子数密度与BU工作流任务单.md` and `2026-05-01_PNJL可选功能盘点与优先级任务单.md`.

## Written / Configured

- Formal plotting entrypoint: `scripts/analysis/pnjl/build_figure4_phase_diagram_assets.py`.
- Result-side figure assets and summary JSON.
- Figure-side PNG/PDF and `plot_manifest.json`.
- Result-side README, audit, and manifest.

## Effective / Usable

- Formal phase reference is production-grade for the stated PNJL phase-reference scope.
- Formal default figure is generated from result-side formal reference data, not old PNGs or paper-side candidates.
- `T-mu_B` axis conversion is explicit and recorded in both result-side assets and figure-side manifest.
- Spinodal data are production data and the plotting script supports an opt-in spinodal layer.

## Not Run / Skipped

- No full refined C2 run for every `xi=-0.5:0.05:0.5` point; C2 is an anchor-grid refined check by design.
- No paper-project integration was run because it is outside this production case scope.

## Known Limitations And Residual Risks

- C2 refined check covers anchor `xi` values rather than a full refined grid.
- Crossover `rho` and derivative diagnostics have larger relative differences at some points; core `T-mu_B` phase-line observables remain stable.
- Spinodal layer can visually crowd the figure, so it is off by default.
- Default plotted subset is a design choice for readability; future figure variants can use the full `xi` grid from the same formal data.
- CEP position in the `T-rho` panel uses a plotting coordinate estimated from nearest boundary densities; visual connector rows use that coordinate for line continuity. The exact CEP source remains the CEP CSV and `T-mu_B` coordinate.
