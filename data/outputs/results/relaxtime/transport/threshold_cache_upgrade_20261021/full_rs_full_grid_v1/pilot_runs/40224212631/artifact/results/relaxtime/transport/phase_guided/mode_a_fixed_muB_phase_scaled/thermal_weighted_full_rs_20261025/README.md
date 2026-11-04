# thermal_weighted_full_rs_20261025

phase-guided transport canonical case for `mode_a_fixed_muB_phase_scaled`.

## Scope
- mode: `mode_a_fixed_muB_phase_scaled`
- case name: `thermal_weighted_full_rs_20261025`
- summary: 固定 muB，沿相变参考温度做 T/T_phase 倍率带，并对每个倍率带连续扫描 xi。
- xi list: `-0.04, -0.003`
- muB list (MeV): `900.0`
- alpha_T list: `1.0`
- compute bulk viscosity (`zeta`): `true`
- propagator xi policy: `match_thermo`
- sigma cache policy: `thermal_weighted`
- thermodynamic nodes: `p_num=96, t_num=24`
- phase anchor policy: `direct_coexistence`
- phase reference root: `/home/runner/work/Julia_RelaxTime/Julia_RelaxTime/data/reference/pnjl/issue130_phase_reference_v2`
- phase reference layer: `accepted`
- phase reference mode: `runtime`
- phase reference source: `candidate`
- phase reference primary layer: `accepted`
- phase reference runtime view: `accepted_primary`
- tau/sigma overrides:
  - tau_p_nodes: `48`
  - tau_angle_nodes: `36`
  - tau_phi_nodes: `36`
  - tau_n_sigma_points: `24`
- channel diagnostics: `true`
- total planned points: `2`

## Key Files
- `sampling_plan.csv`: phase-guided sampling plan
- `phase_guided_transport_scan.csv`: transport scan output CSV
- `effective_config.json`: effective run config snapshot
- `run_manifest.json`: provenance metadata
- canonical figures: `/home/runner/work/Julia_RelaxTime/Julia_RelaxTime/data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/thermal_weighted_full_rs_20261025`

## Interpretation Boundary
- This directory is a user-facing result asset, not an external validation truth set.
- Numerical drift should be governed separately by regression coverage.
- For a directly anchored first-order alpha_T=1 slice, xi=0 is intentionally absent; certified negative/positive near-zero points represent the quark/hadron side limits.
