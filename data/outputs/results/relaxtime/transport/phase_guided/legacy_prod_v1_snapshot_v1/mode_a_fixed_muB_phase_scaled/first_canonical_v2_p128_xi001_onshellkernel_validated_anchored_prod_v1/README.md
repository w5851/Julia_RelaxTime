# first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1 (mode A)

This is the candidate replacement for the superseded RS-shared-energy production case after Issue #130.

## Status

- numerical verdict: `production-grade`
- production registry status: `current_candidate`
- manuscript eligible: `false` until PR review approves promotion

## Physics policies

- transport kernel energy: `onshell_kinematic_energy`, E_kin=sqrt(p^2+M^2)
- distribution energy: `rs_deformed`, E_dist=sqrt(p^2+M^2+xi*(p*cos(theta))^2)
- propagator xi policy: `match_thermo`
- sigma cache policy: `validated_anchored`

## Scope

- mode: `mode_a_fixed_muB_phase_scaled` (fixed mu_B with phase-scaled temperature bands)
- xi grid: `-0.50:0.01:0.50`
- rows: `909`
- source commit: `05be2c05186f8e12baf3097b68f8619e53d19711`
- source candidate case: `p128_xi001_onshellkernel_localdispatch_20260715`

## Evidence

- p104-to-p128 gate: `convergence/parameter_gate/convergence_summary.json`
- full-grid old/new comparison: `convergence/old_vs_new_full_grid_comparison.csv`
- production import summary: `convergence/production_import_summary.json`
- result manifest: `manifest.json`
- figure directory: `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1`
- plot manifest: `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_manifest.json`

The old production directories are retained byte-for-byte and are not modified by this case.

## Generated Figures
- `tau_u`, `tau_d`, `tau_s`, `tau_ubar`, `tau_dbar`, `tau_sbar`, `eta`, `sigma`, `zeta`, `eta_over_s`, `zeta_over_s`, `sigma_over_T` all plot against `xi`
- panel rule: `fixed mu_B panel`
- line rule: `alpha_T lines with explicit T`
- plot manifest: `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_manifest.json`

Generated PNG files:
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB0.0/eta_over_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB0.0/eta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB0.0/sigma_over_T_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB0.0/sigma_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB0.0/tau_d_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB0.0/tau_dbar_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB0.0/tau_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB0.0/tau_sbar_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB0.0/tau_u_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB0.0/tau_ubar_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB0.0/zeta_over_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB0.0/zeta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB450.0/eta_over_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB450.0/eta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB450.0/sigma_over_T_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB450.0/sigma_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB450.0/tau_d_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB450.0/tau_dbar_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB450.0/tau_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB450.0/tau_sbar_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB450.0/tau_u_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB450.0/tau_ubar_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB450.0/zeta_over_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB450.0/zeta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB900.0/eta_over_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB900.0/eta_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB900.0/sigma_over_T_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB900.0/sigma_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB900.0/tau_d_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB900.0/tau_dbar_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB900.0/tau_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB900.0/tau_sbar_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB900.0/tau_u_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB900.0/tau_ubar_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB900.0/zeta_over_s_vs_xi.png`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_panel=muB900.0/zeta_vs_xi.png`
