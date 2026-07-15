# Production Audit: first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1 (mode A)

## Verdict

`production-grade` numerical artifact; registry status remains `current_candidate` pending human review.

## Production Case

- case: `first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1`
- mode: `mode_a_fixed_muB_phase_scaled`
- source commit: `05be2c05186f8e12baf3097b68f8619e53d19711`
- source candidate: `p128_xi001_onshellkernel_localdispatch_20260715`
- workflow: `.github/workflows/relaxtime-phase-guided-transport-production.yml`

## Physics Scope

- `fixed mu_B with phase-scaled temperature bands`
- xi grid `-0.50:0.01:0.50`
- eta, sigma, zeta and their normalized ratios, plus six relaxation times and inverse times
- on-shell E_kin in explicit transport kernels; RS-deformed E_dist only in distributions

## Non-goals

- no change to s, v_rel, t limits, Pauli blocking, or relaxation-time semantics
- no longitudinal/transverse transport tensor decomposition
- this workflow does not output kappa_XY or lambda columns

## Selected Production Parameters

- tau nodes: `128/20/36/24`
- sigma grid: `560`
- sigma cache: `validated_anchored`
- propagator xi policy: `match_thermo`
- bulk viscosity: enabled
- channel diagnostics: enabled

## Convergence And Import Gates

- p104-to-p128 gate verdict: `production-grade`
- unchanged-contract worst relative drift: `0.00016650389`
- tau old/new worst relative drift: `0.00016650389`
- xi=0 transport worst relative drift: `3.1101331e-07`
- source failed rows: `0`
- scan non-finite values: `0`
- diagnostic non-finite values: `0`
- diagnostic negative rate/contribution values: `0`

## Source Action Runs

| shard | run |
| --- | --- |
| modeA_muB0_alpha1.0 / fixed-muB-phase-scaled | [29398949410](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398949410) |
| modeA_muB0_alpha1.1 / fixed-muB-phase-scaled | [29398949754](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398949754) |
| modeA_muB0_alpha1.2 / fixed-muB-phase-scaled | [29398949179](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398949179) |
| modeA_muB450_alpha1.0_xi_neg / fixed-muB-phase-scaled | [29398949133](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398949133) |
| modeA_muB450_alpha1.0_xi_pos / fixed-muB-phase-scaled | [29398949153](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398949153) |
| modeA_muB450_alpha1.1_xi_neg / fixed-muB-phase-scaled | [29398958778](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398958778) |
| modeA_muB450_alpha1.1_xi_pos / fixed-muB-phase-scaled | [29398958810](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398958810) |
| modeA_muB450_alpha1.2_xi_neg / fixed-muB-phase-scaled | [29398958832](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398958832) |
| modeA_muB450_alpha1.2_xi_pos / fixed-muB-phase-scaled | [29398958794](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398958794) |
| modeA_muB900_alpha1.0_xi_neg / fixed-muB-phase-scaled | [29398958714](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398958714) |
| modeA_muB900_alpha1.0_xi_pos / fixed-muB-phase-scaled | [29398966564](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398966564) |
| modeA_muB900_alpha1.1_xi_neg / fixed-muB-phase-scaled | [29398966592](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398966592) |
| modeA_muB900_alpha1.1_xi_pos / fixed-muB-phase-scaled | [29398966543](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398966543) |
| modeA_muB900_alpha1.2_xi_neg / fixed-muB-phase-scaled | [29398966479](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398966479) |
| modeA_muB900_alpha1.2_xi_pos / fixed-muB-phase-scaled | [29398966574](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398966574) |

## Data Outputs

- `phase_guided_transport_scan.csv`
- `sampling_plan.csv`
- `channel_diagnostics.csv`
- `failed_points.csv`
- `effective_config.json`
- `convergence/`

## Figure Outputs

- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_manifest.json`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/**/*.png`

## Known Limitations And Residual Risks

- Action outputs were `diagnostic-only`; this local audit is the explicit production gate.
- Dense-xi local structures still require neighboring-point and channel-diagnostic interpretation.
- Manuscript eligibility remains false until the independent production PR is reviewed.

## Final Local Validation And Command Log

- mode: `mode_a_fixed_muB_phase_scaled`.
- Action collector: `collect_phase_guided_action_artifacts.py --case-name p128_xi001_onshellkernel_localdispatch_20260715 --source-commit 05be2c05186f8e12baf3097b68f8619e53d19711 --expected-count 30 --download`.
- Import gate: `import_phase_guided_energy_semantics_artifacts.py ... --new-case first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1 --validate-only` passed before repository writes.
- Formal import: the same importer without `--validate-only` wrote the new slug only.
- Figure command: `julia --project=. scripts/relaxtime/run_phase_guided_transport_plots.jl --case-dir data/outputs/results/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1 --overwrite`.
- final scan rows: `909`; diagnostic rows: `38178`; converged rows: `909`.
- failed rows: `0`; scan NaN/Inf: `0`; scan negative values: `0`.
- diagnostic NaN/Inf: `0`; diagnostic negative rate/contribution values: `0`.
- figure manifest: `36` PNG SHA256 entries; hash verification: `passed`.
- plot DPI target: `600`; PNG pHYs metadata range: `599.9988` to `599.9988` dpi. The sub-600 display is PNG pixels-per-meter quantization of the 600 dpi target.
- figure-side hygiene: `passed`; only PNG files and `plot_manifest.json` remain after provenance sidecar cleanup.
- result and figure manifests are refreshed after this section is written.
- old-case byte-integrity is checked independently with pre/post tree SHA256 and recorded in the M2 task document and production PR.
