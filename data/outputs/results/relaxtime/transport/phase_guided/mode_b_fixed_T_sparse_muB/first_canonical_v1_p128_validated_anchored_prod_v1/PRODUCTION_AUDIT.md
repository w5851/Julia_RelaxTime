# Production Audit: Phase-Guided Transport mode B p128 validated-anchored

Verdict: `production-grade`

Date: 2026-06-10

## Production Case

- Case slug: `first_canonical_v1_p128_validated_anchored_prod_v1`
- Mode: `mode_b_fixed_T_sparse_muB`
- Result path: `data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1_p128_validated_anchored_prod_v1/`
- Figure path: `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1_p128_validated_anchored_prod_v1/`
- Remote workflow: `.github/workflows/relaxtime-phase-guided-transport-production.yml`
- Source commit: `c3f25dde5c9e62e7a3ebcd8ed50aef2b53377c0f`

## Physics Scope

- Fixed temperature, sparse `mu_B`, continuous `xi` scan.
- Grid: `T = 120, 160, 200 MeV`; `mu_B = 0, 450, 900 MeV`; `xi = -0.50:0.05:0.50`.
- Observables: relaxation times, inverse relaxation times, `eta`, `sigma`, `zeta`, `eta_over_s`, `zeta_over_s`, `sigma_over_T`.
- Bulk viscosity computation: enabled.
- Channel diagnostics: enabled.

## Policies

- `propagator_xi_policy = match_thermo`
- `sigma_cache_policy = validated_anchored`
- `channel_diagnostics = true`

The `validated_anchored` sigma-cache policy is a numerical production policy. It addresses the threshold-subtraction / sigma-cache sensitivity found in the older `first_canonical_v1` artifact while retaining the physical `match_thermo` propagator-xi policy.

## Command Log

The p128 production artifact was generated with GitHub Actions shards using:

- `case_name = first_canonical_v1_match_validated_anchored_plus5_candidate_p128_a20_phi36_ns24_g560`
- `propagator_xi_policy = match_thermo`
- `sigma_cache_policy = validated_anchored`
- `tau_p_nodes = 128`
- `tau_angle_nodes = 20`
- `tau_phi_nodes = 36`
- `tau_n_sigma_points = 24`
- `sigma_grid_n = 560`
- `render_plots = true`

The formal repository case imports the merged artifact under `first_canonical_v1_p128_validated_anchored_prod_v1`. Path metadata in `effective_config.json`, CSV preambles, and plot manifests was rewritten to the repository production paths; numerical CSV rows were not recomputed during import.

Mode-B p128 shard runs:

| shard | run id |
| --- | --- |
| `b_T120_muB0` | `27246617302` |
| `b_T120_muB450` | `27246622499` |
| `b_T120_muB900` | `27246627637` |
| `b_T160_muB0` | `27246634399` |
| `b_T160_muB450` | `27246639390` |
| `b_T160_muB900` | `27246646432` |
| `b_T200_muB0` | `27246651646` |
| `b_T200_muB450` | `27246657527` |
| `b_T200_muB900` | `27246669072` |

## Convergence Matrix

Shared convergence evidence:

- `data/outputs/results/relaxtime/transport/phase_guided/first_canonical_v1_p128_validated_anchored_prod_v1_convergence/p104_vs_p128_convergence_gate.csv`
- `data/outputs/results/relaxtime/transport/phase_guided/first_canonical_v1_p128_validated_anchored_prod_v1_convergence/p104_vs_p128_convergence_gate.summary.json`
- `data/outputs/results/relaxtime/transport/phase_guided/first_canonical_v1_p128_validated_anchored_prod_v1_convergence/p104_vs_p128_convergence_gate.md`

Adjacent high-precision comparison:

| parameter set | role |
| --- | --- |
| `104/18/30/20/440` | adjacent reference |
| `128/20/36/24/560` | selected production |

Summary:

| metric | value |
| --- | ---: |
| matched points | `378` |
| overall max relative difference | `0.009426486578504083` |
| core max relative difference | `0.006176952435828459` |
| anti-particle max relative difference | `0.009426486578504083` |
| count above 1% | `0` |
| count above 0.5% | `12` |

Worst p104-to-p128 point:

- `mode_a`, `muB=450 MeV`, `alpha_T=1.1`, `xi=0.35`, `tau_ubar/tau_dbar`, relative difference `0.9426486578504083%`

Largest listed mode-B watch window:

- `mode_b`, `T=120 MeV`, `muB=900 MeV`, `xi=-0.15`, max relative difference `0.154604%`

## Data Outputs

- `phase_guided_transport_scan.csv`: `189` data rows
- `sampling_plan.csv`: `189` data rows
- `channel_diagnostics.csv`: `7938` data rows
- `failed_points.csv`: `0` failure rows
- `effective_config.json`
- `README.md`
- `PRODUCTION_AUDIT.md`
- `manifest.json`

Channel diagnostics audit:

- non-finite diagnostic rows: `0`
- negative rate/contribution rows: `0`

## Figure Outputs

- Figure directory: `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1_p128_validated_anchored_prod_v1/`
- PNG files: `24`
- `plot_manifest.json`: present
- smallest PNG size: `146712` bytes

## Old Result Supersession

This case supersedes:

- `data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1/`

The old result is retained as a historical artifact. It lacks explicit high-precision tau/sigma integration parameters, `validated_anchored` sigma-cache treatment, channel diagnostics, and a full-grid convergence gate. It should not be used for production-grade interpretation of local `xi` structures.

## Validation Commands And Results

Repository import checks are recorded in the PR validation log. The production-grade verdict is based on the action artifact audit plus the p104-to-p128 convergence gate summarized above.

## Known Limitations And Residual Risks

- The 1% gate is an adjacent high-precision convergence criterion, not a mathematical proof of exact convergence.
- Local structures that remain after the `validated_anchored` correction should be interpreted through the physical `match_thermo` propagator and denominator-chain diagnostics, not as evidence of the old threshold-subtraction cache artifact.
- The older `first_canonical_v1` result is not deleted in this update to preserve provenance and comparison ability.
