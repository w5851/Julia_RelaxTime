# first_canonical_v1

phase-guided transport canonical case for `mode_b_fixed_T_sparse_muB`.

## Scope
- mode: `mode_b_fixed_T_sparse_muB`
- case name: `first_canonical_v1`
- summary: 固定温度、离散 muB、连续 xi 的相图邻域稀疏扫描。
- xi list: `-0.5, -0.2, 0.0, 0.2`
- muB list (MeV): `0.0, 450.0, 900.0`
- fixed T list (MeV): `120.0, 160.0, 200.0`
- compute bulk viscosity (`zeta`): `true`
- total planned points: `36`

## Key Files
- `sampling_plan.csv`: phase-guided sampling plan
- `phase_guided_transport_scan.csv`: transport scan output CSV
- `effective_config.json`: effective run config snapshot
- `run_manifest.json`: provenance metadata
- canonical figures: `D:/Desktop/Julia_RelaxTime/data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1`

## Interpretation Boundary
- This directory is a user-facing result asset, not an external validation truth set.
- Numerical drift should be governed separately by regression coverage.

## Generated Figures
- `tau_u`, `tau_d`, `tau_s`, `eta`, `sigma`, `zeta`, `eta_over_s`, `sigma_over_T` all plot against `xi`
- panel rule: `fixed T panel`
- line rule: `mu_B lines`
- plot manifest: `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_manifest.json`

Generated PNG files:
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T120.0/eta_over_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T120.0/eta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T120.0/sigma_over_T_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T120.0/sigma_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T120.0/tau_d_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T120.0/tau_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T120.0/tau_u_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T120.0/zeta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T160.0/eta_over_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T160.0/eta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T160.0/sigma_over_T_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T160.0/sigma_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T160.0/tau_d_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T160.0/tau_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T160.0/tau_u_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T160.0/zeta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T200.0/eta_over_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T200.0/eta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T200.0/sigma_over_T_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T200.0/sigma_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T200.0/tau_d_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T200.0/tau_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T200.0/tau_u_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/plot_panel=T200.0/zeta_vs_xi.png`
