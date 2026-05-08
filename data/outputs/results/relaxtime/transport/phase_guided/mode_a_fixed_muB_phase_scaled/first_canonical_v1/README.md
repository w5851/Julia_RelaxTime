# first_canonical_v1

phase-guided transport canonical case for `mode_a_fixed_muB_phase_scaled`.

## Scope
- mode: `mode_a_fixed_muB_phase_scaled`
- case name: `first_canonical_v1`
- summary: 固定 muB，沿相变参考温度做 T/T_phase 倍率带，并对每个倍率带连续扫描 xi。
- xi list: `-0.5, -0.2, 0.0, 0.2`
- muB list (MeV): `260.0, 400.0`
- alpha_T list: `1.0, 1.1`
- compute bulk viscosity (`zeta`): `true`
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

Generated PNG files:
- `figures/scan_group=muB260.0_alpha1.0/eta_over_s_vs_xi.png`
- `figures/scan_group=muB260.0_alpha1.0/eta_vs_xi.png`
- `figures/scan_group=muB260.0_alpha1.0/multi_y_tau_u_tau_d_tau_s_vs_xi.png`
- `figures/scan_group=muB260.0_alpha1.0/sigma_over_T_vs_xi.png`
- `figures/scan_group=muB260.0_alpha1.0/sigma_vs_xi.png`
- `figures/scan_group=muB260.0_alpha1.0/zeta_vs_xi.png`
- `figures/scan_group=muB260.0_alpha1.1/eta_over_s_vs_xi.png`
- `figures/scan_group=muB260.0_alpha1.1/eta_vs_xi.png`
- `figures/scan_group=muB260.0_alpha1.1/multi_y_tau_u_tau_d_tau_s_vs_xi.png`
- `figures/scan_group=muB260.0_alpha1.1/sigma_over_T_vs_xi.png`
- `figures/scan_group=muB260.0_alpha1.1/sigma_vs_xi.png`
- `figures/scan_group=muB260.0_alpha1.1/zeta_vs_xi.png`
- `figures/scan_group=muB400.0_alpha1.0/eta_over_s_vs_xi.png`
- `figures/scan_group=muB400.0_alpha1.0/eta_vs_xi.png`
- `figures/scan_group=muB400.0_alpha1.0/multi_y_tau_u_tau_d_tau_s_vs_xi.png`
- `figures/scan_group=muB400.0_alpha1.0/sigma_over_T_vs_xi.png`
- `figures/scan_group=muB400.0_alpha1.0/sigma_vs_xi.png`
- `figures/scan_group=muB400.0_alpha1.0/zeta_vs_xi.png`
- `figures/scan_group=muB400.0_alpha1.1/eta_over_s_vs_xi.png`
- `figures/scan_group=muB400.0_alpha1.1/eta_vs_xi.png`
- `figures/scan_group=muB400.0_alpha1.1/multi_y_tau_u_tau_d_tau_s_vs_xi.png`
- `figures/scan_group=muB400.0_alpha1.1/sigma_over_T_vs_xi.png`
- `figures/scan_group=muB400.0_alpha1.1/sigma_vs_xi.png`
- `figures/scan_group=muB400.0_alpha1.1/zeta_vs_xi.png`
