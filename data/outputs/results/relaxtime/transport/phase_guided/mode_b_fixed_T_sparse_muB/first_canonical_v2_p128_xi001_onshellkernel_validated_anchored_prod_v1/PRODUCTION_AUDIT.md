# Production Audit: first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1 (mode B)

## Verdict

`production-grade` numerical artifact; registry status remains `current_candidate` pending human review.

## Production Case

- case: `first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1`
- mode: `mode_b_fixed_T_sparse_muB`
- source commit: `05be2c05186f8e12baf3097b68f8619e53d19711`
- source candidate: `p128_xi001_onshellkernel_localdispatch_20260715`
- workflow: `.github/workflows/relaxtime-phase-guided-transport-production.yml`

## Physics Scope

- `fixed temperature with sparse mu_B values`
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
- unchanged-contract worst relative drift: `8.2279794e-06`
- tau old/new worst relative drift: `8.2279794e-06`
- xi=0 transport worst relative drift: `7.9976584e-07`
- source failed rows: `0`
- scan non-finite values: `0`
- diagnostic non-finite values: `0`
- diagnostic negative rate/contribution values: `0`

## Source Action Runs

| shard | run |
| --- | --- |
| modeB_T120_muB0 / fixed-T-sparse-muB | [29398990593](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398990593) |
| modeB_T120_muB450_xi_neg / fixed-T-sparse-muB | [29398990688](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398990688) |
| modeB_T120_muB450_xi_pos / fixed-T-sparse-muB | [29398991616](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398991616) |
| modeB_T120_muB900_xi_neg / fixed-T-sparse-muB | [29398997416](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398997416) |
| modeB_T120_muB900_xi_pos / fixed-T-sparse-muB | [29399027156](https://github.com/w5851/Julia_RelaxTime/actions/runs/29399027156) |
| modeB_T160_muB0 / fixed-T-sparse-muB | [29398990586](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398990586) |
| modeB_T160_muB450_xi_neg / fixed-T-sparse-muB | [29399030206](https://github.com/w5851/Julia_RelaxTime/actions/runs/29399030206) |
| modeB_T160_muB450_xi_pos / fixed-T-sparse-muB | [29398997376](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398997376) |
| modeB_T160_muB900_xi_neg / fixed-T-sparse-muB | [29399033033](https://github.com/w5851/Julia_RelaxTime/actions/runs/29399033033) |
| modeB_T160_muB900_xi_pos / fixed-T-sparse-muB | [29399004218](https://github.com/w5851/Julia_RelaxTime/actions/runs/29399004218) |
| modeB_T200_muB0 / fixed-T-sparse-muB | [29398990633](https://github.com/w5851/Julia_RelaxTime/actions/runs/29398990633) |
| modeB_T200_muB450_xi_neg / fixed-T-sparse-muB | [29399005076](https://github.com/w5851/Julia_RelaxTime/actions/runs/29399005076) |
| modeB_T200_muB450_xi_pos / fixed-T-sparse-muB | [29399006769](https://github.com/w5851/Julia_RelaxTime/actions/runs/29399006769) |
| modeB_T200_muB900_xi_neg / fixed-T-sparse-muB | [29399004410](https://github.com/w5851/Julia_RelaxTime/actions/runs/29399004410) |
| modeB_T200_muB900_xi_pos / fixed-T-sparse-muB | [29399005317](https://github.com/w5851/Julia_RelaxTime/actions/runs/29399005317) |

## Data Outputs

- `phase_guided_transport_scan.csv`
- `sampling_plan.csv`
- `channel_diagnostics.csv`
- `failed_points.csv`
- `effective_config.json`
- `convergence/`

## Figure Outputs

- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/plot_manifest.json`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/**/*.png`

## Known Limitations And Residual Risks

- Action outputs were `diagnostic-only`; this local audit is the explicit production gate.
- Dense-xi local structures still require neighboring-point and channel-diagnostic interpretation.
- Manuscript eligibility remains false until the independent production PR is reviewed.

## Final Local Validation And Command Log

- mode: `mode_b_fixed_T_sparse_muB`.
- Action collector: `collect_phase_guided_action_artifacts.py --case-name p128_xi001_onshellkernel_localdispatch_20260715 --source-commit 05be2c05186f8e12baf3097b68f8619e53d19711 --expected-count 30 --download`.
- Import gate: `import_phase_guided_energy_semantics_artifacts.py ... --new-case first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1 --validate-only` passed before repository writes.
- Formal import: the same importer without `--validate-only` wrote the new slug only.
- Figure command: `julia --project=. scripts/relaxtime/run_phase_guided_transport_plots.jl --case-dir data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1 --overwrite`.
- final scan rows: `909`; diagnostic rows: `38178`; converged rows: `909`.
- failed rows: `0`; scan NaN/Inf: `0`; scan negative values: `0`.
- diagnostic NaN/Inf: `0`; diagnostic negative rate/contribution values: `0`.
- figure manifest: `36` PNG SHA256 entries; hash verification: `passed`.
- plot DPI target: `600`; PNG pHYs metadata range: `599.9988` to `599.9988` dpi. The sub-600 display is PNG pixels-per-meter quantization of the 600 dpi target.
- figure-side hygiene: `passed`; only PNG files and `plot_manifest.json` remain after provenance sidecar cleanup.
- result and figure manifests are refreshed after this section is written.
- old-case byte-integrity is checked independently with pre/post tree SHA256 and recorded in the M2 task document and production PR.
