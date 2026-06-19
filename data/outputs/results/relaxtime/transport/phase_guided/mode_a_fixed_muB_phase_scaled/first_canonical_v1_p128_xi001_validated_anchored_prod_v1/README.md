# first_canonical_v1_p128_xi001_validated_anchored_prod_v1

Phase-guided transport production-grade case for `mode_a_fixed_muB_phase_scaled`.

## Verdict

`production-grade`

This result imports the p128 high-xi-resolution GitHub Actions artifacts after local merge gates passed. The numerical integration policy is unchanged from `first_canonical_v1_p128_validated_anchored_prod_v1`; the new contribution is denser `xi = -0.50:0.01:0.50` sampling.

## Scope

- mode: `mode_a_fixed_muB_phase_scaled`
- case name: `first_canonical_v1_p128_xi001_validated_anchored_prod_v1`
- summary: fixed mu_B, phase-reference T/T_phase bands, continuous xi scan
- grid: `mu_B = [0.0, 450.0, 900.0] MeV, alpha_T = [1.0, 1.1, 1.2], xi = -0.50:0.01:0.50`
- planned and completed points: `909`
- compute bulk viscosity (`zeta`): `true`

## Numerical Policy

- `propagator_xi_policy = match_thermo`
- `sigma_cache_policy = validated_anchored`
- `tau_p_nodes = 128`
- `tau_angle_nodes = 20`
- `tau_phi_nodes = 36`
- `tau_n_sigma_points = 24`
- `sigma_grid_n = 560`
- `channel_diagnostics = true`

## Key Files

- `phase_guided_transport_scan.csv`: imported transport scan output
- `sampling_plan.csv`: full phase-guided sampling plan
- `channel_diagnostics.csv`: per-channel scattering-rate diagnostics
- `failed_points.csv`: failure ledger; contains `0` data rows
- `effective_config.json`: repository-imported effective configuration snapshot
- `manifest.json`: production manifest with source runs, hashes, and gate evidence
- `PRODUCTION_AUDIT.md`: human-readable production audit

Shared high-xi import evidence is stored in:

- `data/outputs/results/relaxtime/transport/phase_guided/first_canonical_v1_p128_xi001_validated_anchored_prod_v1_convergence/`

Canonical figures are stored in:

- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/`

## Gate Summary

- source run head SHA: `700a0845abd09013eb39ff8f8d2993d5574476d3`
- source action artifacts: `15` result artifacts
- failed point rows across source artifacts: `0`
- chunk overlap duplicate groups: `6`
- chunk overlap worst relative difference: `4.014249626836331e-06`
- anchor comparison to `first_canonical_v1_p128_validated_anchored_prod_v1` effective worst relative difference: `8.453069929263162e-05`
- anchor comparison raw near-zero worst relative difference: `1.125`
- inherited adjacent integration convergence: `first_canonical_v1_p128_validated_anchored_prod_v1_convergence/p104_vs_p128_convergence_gate.summary.json`

## Interpretation Boundary

This case supersedes neither `first_canonical_v1_p128_validated_anchored_prod_v1` nor its convergence evidence. It is the same p128 validated-anchored numerical policy sampled on a denser `xi` grid. Claims that depend on newly resolved local structures should still check neighboring points and channel diagnostics.

## Generated Figures
- `tau_u`, `tau_d`, `tau_s`, `eta`, `sigma`, `zeta`, `eta_over_s`, `sigma_over_T` all plot against `xi`
- panel rule: `fixed mu_B panel`
- line rule: `alpha_T lines with explicit T`
- plot manifest: `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_manifest.json`

Generated PNG files:
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB0.0/eta_over_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB0.0/eta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB0.0/sigma_over_T_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB0.0/sigma_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB0.0/tau_d_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB0.0/tau_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB0.0/tau_u_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB0.0/zeta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB450.0/eta_over_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB450.0/eta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB450.0/sigma_over_T_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB450.0/sigma_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB450.0/tau_d_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB450.0/tau_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB450.0/tau_u_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB450.0/zeta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB900.0/eta_over_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB900.0/eta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB900.0/sigma_over_T_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB900.0/sigma_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB900.0/tau_d_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB900.0/tau_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB900.0/tau_u_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/plot_panel=muB900.0/zeta_vs_xi.png`
