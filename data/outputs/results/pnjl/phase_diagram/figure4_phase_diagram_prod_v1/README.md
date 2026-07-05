# figure4_phase_diagram_prod_v1

Formal PNJL phase-reference and Figure 4 asset package for anisotropic phase-diagram analysis.

## Verdict

- phase-reference verdict: `production-grade`
- convergence evidence: `convergence/P4-C_CONVERGENCE_SUMMARY.md`
- formal reference tag: `figure4_phase_diagram_prod_v1_c1_p24t8`
- formal reference root: `reference/`

The formal reference covers `xi=-0.5:0.05:0.5` and includes first-order boundary, CEP, crossover, and spinodal data. The default figure uses a readable subset `xi=-0.5,-0.25,0.0,0.25,0.5`; the full-grid data remain available for future figure variants.

## Key Files

- `phase_reference_source_manifest.json`: promoted source manifest for the production reference.
- `reference/boundary_figure4_phase_diagram_prod_v1_c1_p24t8.csv`: first-order boundary in quark chemical potential `muq`.
- `reference/cep_figure4_phase_diagram_prod_v1_c1_p24t8.csv`: CEP table with direct `muB_CEP_MeV`.
- `reference/crossover_figure4_phase_diagram_prod_v1_c1_p24t8.csv`: crossover line in quark chemical potential `mu_MeV`.
- `reference/spinodals_figure4_phase_diagram_prod_v1_c1_p24t8.csv`: spinodal branches in quark chemical potential.
- `figure_assets/figure4_phase_lines_TmuB.csv`: normalized plot asset for the `T-mu_B` panel.
- `figure_assets/figure4_phase_lines_Trho.csv`: normalized plot asset for the `T-rho` panel.
- `figure_assets/figure4_phase_plot_inputs_summary.json`: plot-input summary and row counts.
- `PRODUCTION_AUDIT.md`: production audit and validation record.
- `manifest.json`: result-side file index and hashes.

Figure-side files are written under `data/outputs/figures/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/`:

- `figure4_phase_diagram_TmuB_Trho.png`
- `figure4_phase_diagram_TmuB_Trho.pdf`
- `plot_manifest.json`

## Plot Command

```powershell
python scripts/analysis/pnjl/build_figure4_phase_diagram_assets.py
```

Default behavior:

- source: this package's `reference/` root and tag `figure4_phase_diagram_prod_v1_c1_p24t8`
- plotted subset: `xi=-0.5,-0.25,0.0,0.25,0.5`
- spinodal layer: off in the default PNG/PDF
- outputs: result-side `figure_assets/`, figure-side PNG/PDF and `plot_manifest.json`

To generate a spinodal-layer variant, pass `--include-spinodal` and choose a separate figure directory if the default figure should remain unchanged.

## Axis And Unit Notes

- `T-mu_B` assets explicitly convert source quark chemical potentials with `muB_MeV = 3 * muq_MeV`.
- CEP uses `muB_CEP_MeV` directly when present.
- `T-rho` uses `rho/rho0`.
- CEP `rho` in the `T-rho` asset is an estimated plotting coordinate: the mean of boundary-side densities at the nearest available boundary temperature for the same `xi`.
- `T-rho` contains plotting-only CEP connector rows so that the coexistence branches and crossover curve visually pass through CEP. The formal source reference CSV files are unchanged.

## Consumption Boundary

This package is a Julia_RelaxTime formal production artifact. It does not modify the paper project, captions, claim tables, LaTeX, or manuscript narrative.
