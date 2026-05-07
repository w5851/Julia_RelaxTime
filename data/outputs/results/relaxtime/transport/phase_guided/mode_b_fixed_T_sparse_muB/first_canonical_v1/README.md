# first_canonical_v1

phase-guided transport canonical case for `mode_b_fixed_T_sparse_muB`.

## Scope
- mode: `mode_b_fixed_T_sparse_muB`
- case name: `first_canonical_v1`
- summary: 固定温度、离散 muB、连续 xi 的相图邻域稀疏扫描。
- xi list: `-0.5, -0.2, 0.0, 0.2`
- muB list (MeV): `0.0, 260.0`
- fixed T list (MeV): `120.0, 138.0`
- total planned points: `16`

## Key Files
- `sampling_plan.csv`: phase-guided sampling plan
- `phase_guided_transport_scan.csv`: transport scan output CSV
- `effective_config.json`: effective run config snapshot
- `run_manifest.json`: provenance metadata
- `figures/`: reserved for plot-review / canonical figures

## Interpretation Boundary
- This directory is a user-facing result asset, not an external validation truth set.
- Numerical drift should be governed separately by regression coverage.

## Generated Figures
- tau flavor-resolved same-panel figures by `scan_group`
- individual `eta`, `sigma`, `zeta`, `eta_over_s`, `sigma_over_T` figures by `scan_group`
- plot manifest: `figures/plot_manifest.json`
- current manifest count: `20` generated PNG files
- `zeta` / `zeta_over_s` are present in the CSV contract, but this canonical case currently yields `NaN` across all sampled points, so no `zeta` PNG is emitted by the plotting pass

Generated PNG files:
- `figures/scan_group=T120.0_muB0.0/eta_over_s_vs_xi.png`
- `figures/scan_group=T120.0_muB0.0/eta_vs_xi.png`
- `figures/scan_group=T120.0_muB0.0/multi_y_tau_u_tau_d_tau_s_vs_xi.png`
- `figures/scan_group=T120.0_muB0.0/sigma_over_T_vs_xi.png`
- `figures/scan_group=T120.0_muB0.0/sigma_vs_xi.png`
- `figures/scan_group=T120.0_muB260.0/eta_over_s_vs_xi.png`
- `figures/scan_group=T120.0_muB260.0/eta_vs_xi.png`
- `figures/scan_group=T120.0_muB260.0/multi_y_tau_u_tau_d_tau_s_vs_xi.png`
- `figures/scan_group=T120.0_muB260.0/sigma_over_T_vs_xi.png`
- `figures/scan_group=T120.0_muB260.0/sigma_vs_xi.png`
- `figures/scan_group=T138.0_muB0.0/eta_over_s_vs_xi.png`
- `figures/scan_group=T138.0_muB0.0/eta_vs_xi.png`
- `figures/scan_group=T138.0_muB0.0/multi_y_tau_u_tau_d_tau_s_vs_xi.png`
- `figures/scan_group=T138.0_muB0.0/sigma_over_T_vs_xi.png`
- `figures/scan_group=T138.0_muB0.0/sigma_vs_xi.png`
- `figures/scan_group=T138.0_muB260.0/eta_over_s_vs_xi.png`
- `figures/scan_group=T138.0_muB260.0/eta_vs_xi.png`
- `figures/scan_group=T138.0_muB260.0/multi_y_tau_u_tau_d_tau_s_vs_xi.png`
- `figures/scan_group=T138.0_muB260.0/sigma_over_T_vs_xi.png`
- `figures/scan_group=T138.0_muB260.0/sigma_vs_xi.png`
