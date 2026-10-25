# thermal_weighted_full_rs_mode_b_20261025

phase-guided transport canonical case for `mode_b_fixed_T_sparse_muB`.

## Scope
- mode: `mode_b_fixed_T_sparse_muB`
- case name: `thermal_weighted_full_rs_mode_b_20261025`
- summary: 固定温度、离散 muB、连续 xi 的相图邻域稀疏扫描。
- xi list: `-0.5, -0.13, -0.12`
- muB list (MeV): `900.0`
- fixed T list (MeV): `120.0`
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
- total planned points: `3`

## Key Files
- `sampling_plan.csv`: phase-guided sampling plan
- `phase_guided_transport_scan.csv`: transport scan output CSV
- `effective_config.json`: effective run config snapshot
- `run_manifest.json`: provenance metadata
- canonical figures: `/home/runner/work/Julia_RelaxTime/Julia_RelaxTime/data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/thermal_weighted_full_rs_mode_b_20261025`

## Interpretation Boundary
- This directory is a user-facing result asset, not an external validation truth set.
- Numerical drift should be governed separately by regression coverage.
- For a directly anchored first-order alpha_T=1 slice, xi=0 is intentionally absent; certified negative/positive near-zero points represent the quark/hadron side limits.
