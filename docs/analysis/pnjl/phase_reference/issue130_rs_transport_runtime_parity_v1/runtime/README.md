# runtime

phase-guided transport canonical case for `mode_a_fixed_muB_phase_scaled`.

## Scope
- mode: `mode_a_fixed_muB_phase_scaled`
- case name: `issue130_rs_parity_v1`
- summary: 固定 muB，沿相变参考温度做 T/T_phase 倍率带，并对每个倍率带连续扫描 xi。
- xi list: `0.0`
- muB list (MeV): `150.0`
- alpha_T list: `1.0`
- compute bulk viscosity (`zeta`): `false`
- propagator xi policy: `isotropic`
- sigma cache policy: `validated_anchored`
- thermodynamic nodes: `p_num=12, t_num=6`
- phase anchor policy: `direct_coexistence`
- phase reference root: `D:\Temp\julia_relaxtime_issue130_shared\data\reference\pnjl\issue130_phase_reference_v1`
- phase reference layer: `strict`
- phase reference mode: `runtime`
- phase reference source: `candidate`
- phase reference runtime view: `certified_candidate_with_legacy_fallback`
- tau/sigma overrides:
  - tau_p_nodes: `6`
  - tau_angle_nodes: `2`
  - tau_phi_nodes: `2`
  - tau_n_sigma_points: `4`
  - sigma_grid_n: `12`
- channel diagnostics: `false`
- total planned points: `1`

## Key Files
- `sampling_plan.csv`: phase-guided sampling plan
- `phase_guided_transport_scan.csv`: transport scan output CSV
- `effective_config.json`: effective run config snapshot
- `run_manifest.json`: provenance metadata
- canonical figures: `D:/Temp/julia_relaxtime_issue130_shared/data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/issue130_rs_parity_v1`

## Interpretation Boundary
- This directory is a user-facing result asset, not an external validation truth set.
- Numerical drift should be governed separately by regression coverage.
- For a directly anchored first-order alpha_T=1 slice, xi=0 is intentionally absent; certified negative/positive near-zero points represent the quark/hadron side limits.
