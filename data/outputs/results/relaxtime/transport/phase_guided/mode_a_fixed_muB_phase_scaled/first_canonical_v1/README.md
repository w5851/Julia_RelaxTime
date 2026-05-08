# first_canonical_v1

phase-guided transport canonical case for `mode_a_fixed_muB_phase_scaled`.

## Scope
- mode: `mode_a_fixed_muB_phase_scaled`
- case name: `first_canonical_v1`
- summary: 固定 muB，沿相变参考温度做 T/T_phase 倍率带，并对每个倍率带连续扫描 xi。
- xi list: `-0.5, -0.2, 0.0, 0.2`
- muB list (MeV): `0.0, 450.0, 900.0`
- alpha_T list: `1.0, 1.1, 1.2`
- compute bulk viscosity (`zeta`): `true`
- total planned points: `36`

## Key Files
- `sampling_plan.csv`: phase-guided sampling plan
- `phase_guided_transport_scan.csv`: transport scan output CSV
- `effective_config.json`: effective run config snapshot
- `run_manifest.json`: provenance metadata
- canonical figures: `D:/Desktop/Julia_RelaxTime/data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1`

## Interpretation Boundary
- This directory is a user-facing result asset, not an external validation truth set.
- Numerical drift should be governed separately by regression coverage.

## Generated Figures
- `tau_u`, `tau_d`, `tau_s`, `eta`, `sigma`, `zeta`, `eta_over_s`, `sigma_over_T` all plot against `xi`
- panel rule: `fixed mu_B panel`
- line rule: `alpha_T lines`
- plot manifest: `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_manifest.json`

Generated PNG files:
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB0.0/eta_over_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB0.0/eta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB0.0/sigma_over_T_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB0.0/sigma_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB0.0/tau_d_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB0.0/tau_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB0.0/tau_u_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB0.0/zeta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB450.0/eta_over_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB450.0/eta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB450.0/sigma_over_T_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB450.0/sigma_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB450.0/tau_d_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB450.0/tau_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB450.0/tau_u_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB450.0/zeta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB900.0/eta_over_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB900.0/eta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB900.0/sigma_over_T_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB900.0/sigma_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB900.0/tau_d_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB900.0/tau_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB900.0/tau_u_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1/plot_panel=muB900.0/zeta_vs_xi.png`
