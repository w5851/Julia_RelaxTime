# Production Audit: trho_asymmetric K+/pi+ Meson Density Scan v1

Verdict: `blocked`

Date: 2026-06-06

## Production Case

- Case slug: `trho_asymmetric_kplus_piplus_scan_v1`
- Result path: `data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/`
- Figure path: `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/`
- Remote workflow: `.github/workflows/relaxtime-meson-density-production.yml`
- Remote run: `27064914082`
- Source commit: `05d181730bff5d7498918e0c1dceddaf65b5448f`

## Physics Scope

- Equilibrium source: `FixedAsymmetricRho`
- Constraint targets:
  - `sum(rho_i)/(3 rho0) = rho_target`
  - `rho_u / rho_d = 0.876`
  - `rho_s = 0 fm^-3`
- Charged profile: `asymmetric_kplus_over_piplus_signed`
- Meson channels: `pi+`, `K+`
- Density regimes:
  - `stable`
  - `strict_bw_stage1`
  - `phase_shift_current`
  - `phase_shift_gbu_reference`
- Phase policy:
  - `real_axis_mode = pv_b0_eta0`
  - `phase_convention = arg_inverse_propagator`
  - `phase_display = fold_0_pi`
  - `density_policy = strict_normal_domain`
  - `noanom_policy = none`

## Non-Goals

- No new solver mode.
- No change to `FixedAsymmetricRho` constraints.
- No numerical regression baseline update.
- No production-grade data generation after the convergence gate failed.
- No switch to `K- / pi-`, narrowed grid, or diagnostic density policy without a follow-up scope decision.

## Command Log

Remote `convergence_low` command is recorded in:

`convergence/convergence_low/run.command.txt`

The run used:

- `T = 120:10:220 MeV`
- `rho/rho0 = 0.05:0.05:1.00`
- `stable_q_nodes = 96`
- `q_nodes = 12`
- `omega_nodes = 12`
- `omega_min = 0.05 fm^-1`
- `omega_max = 4.0 fm^-1`

## Convergence Matrix

Only `convergence_low` was retained because it already triggered the formal
hard gate. Running `mid` and `high` would refine the same invalid-domain
failure, not establish production-grade convergence.

| item | value |
| --- | ---: |
| rows | 880 |
| path points | 220 |
| all-regime-ok points | 13 |
| rows with `ok` status | 245 |
| rows with `failed` status | 221 |
| rows with `unsafe_bose_domain` status | 414 |
| stable failed rows | 14 |
| strict BW failed rows | 207 |
| phase-shift unsafe rows | 414 |

By regime:

| regime | ok | failed | unsafe_bose_domain |
| --- | ---: | ---: | ---: |
| stable | 206 | 14 | 0 |
| strict_bw_stage1 | 13 | 207 | 0 |
| phase_shift_current | 13 | 0 | 207 |
| phase_shift_gbu_reference | 13 | 0 | 207 |

Constraint diagnostics were exact at the printed precision:

- `constraint_residual_norm`: min `0.0`, max `0.0`
- `rho_u_over_rho_d`: min `0.876`, max `0.876`
- `rho_s_fm3`: min `0.0`, max `0.0`

## Blocker

The selected `K+ / pi+` asymmetric profile produces positive `mu_K`. In the
attempted grid, this creates Bose-domain violations:

- stable density can fail when `mass <= mu_K`;
- strict BW rejects points where `omega_min <= mu_K`;
- phase-shift regimes under `strict_normal_domain` reject support with
  `omega <= mu_K`.

This is a domain-definition problem, not evidence that the quadrature node
count is insufficient.

Follow-up convergence attempts on commit `4c6c1542` applied the explicit
phase-shift-only `density_policy=x_min_cut` policy.

| attempt | remote runs | result |
| --- | --- | --- |
| `bose_x_min=1e-6` | `27112756019`, `27112924358`, `27112924343` | BU/GBU phase-shift rows were runnable, but `mid -> high` K-related quantities still showed order-50% maximum relative differences. |
| `bose_x_min=1e-2` | `27113471270`, `27113471257`, `27113471249`, diagnostic `27113589340` | Cutoff sensitivity improved, but `mid -> high` K-related maximum relative differences remained about 11%--13%; high/custom comparisons exposed upstream `FixedAsymmetricRho` branch changes at the same `(T, rho_target)`. |

The current blocker is therefore upstream density-constrained equilibrium
branch stability in the combined scan path. In particular, some same-grid
points satisfy the printed constraint residual but land on different
`FixedAsymmetricRho` branches across resolution profiles. Production must
remain blocked until the scan path uses a deterministic and documented
continuity/branch policy and the convergence gate passes again.

## Selected Production Parameters

None. No production run was launched because the convergence gate did not reach
`production-grade`.

## Data Outputs

Retained convergence evidence:

- `convergence/convergence_low/combined_meson_density_scan.csv`
- `convergence/convergence_low/README.md`
- `convergence/convergence_low/remote_run_manifest.json`
- `convergence/convergence_low/run.command.txt`
- `convergence/convergence_low/run.stdout.log`
- `convergence/convergence_low/run.stderr.log`
- `convergence/convergence_summary.json`
- `manifest.json`

## Figure Outputs

Retained convergence figure evidence:

- `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/convergence/convergence_low/combined_meson_density_scan.svg`
- `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/convergence/convergence_low/plot_manifest.json`

No production figure was generated.

## Validation Commands And Results

Completed before this audit:

- GitHub Actions remote workflow run `27064914082`: success
- Remote script exit code: `0`

Pending after this audit:

- repository governance checks after committing the blocked-audit artifacts

## Known Limitations And Residual Risks

- `K+ / pi+` under `FixedAsymmetricRho` may require an explicit treatment of
  positive meson chemical potential domains before formal production.
- `trho_asymmetric` formal production depends on upstream
  `FixedAsymmetricRho` branch stability across convergence profiles; this was
  not satisfied in the `bose_x_min` follow-up attempts.
- BW remains a comparison regime only and must not be used alone as a formal
  conclusion.
- A diagnostic policy such as `excitation_only_E_gt_mu` could produce numeric
  tables, but that would be diagnostic-only unless separately justified.
- A channel switch to `K- / pi-` or a narrowed grid is a material scope change
  and requires a follow-up task decision.
