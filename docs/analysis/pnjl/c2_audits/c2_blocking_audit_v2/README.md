# C2 blocking audit v2

这是针对当前 C2 Actions 的 solver-free blocking audit。脚本只读取既有 C1/C2 artifact 和 comparator 输出，不调用 Julia equilibrium solver，不修改正式 reference/result。

## 输入

- C1 run: `31762201725`; tag `issue130_stagec_density_v2_c1_integral_tight_20260813`
- C2 run: `31862752226`; tag `issue130_stagec_density_v2_c2_grid_tight_20260813`
- calculation SHA: `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- postprocess SHA: `fd359e792a89beb5ab12349bba761dc58ee16761`
- comparison policy: fixed public xi grid and physical-key matching from `compare_phase_reference_convergence.py`
- excluded run: `31862709047` (plan failed before numerical jobs; no evidence used)

## 结论

- comparator primary verdict: `classification_regression`
- C2 artifact is accepted as `diagnostic-only convergence candidate`; it is not a phase-reference candidate.
- All retained boundary Maxwell rows are finite/converged, so `rho_geometry_not_converged` is an outer cross-layer certificate failure, not proof that Maxwell bisection failed.
- C2 has `3` CEP bracket(s) above the hard `0.1 MeV` gate in the complete retained CEP table.
- Crossover remains a comparatively stable subsystem under the stored C1/C2 comparison, but this does not close the phase-reference gate.
- At the nine classification-regression anchors, C1 has a boundary row while C2 has no same-temperature boundary row; the audit records this as unavailable geometry, not zero drift.

## Axis and failure distribution

| axis | total | unconverged | fraction |
|---|---:|---:|---:|
| rho | 18511 | 4716 | 25.48% |
| temperature | 4208 | 680 | 16.16% |
| xi | 72 | 28 | 38.89% |

Top failure reasons are stored in `tables/c2_grid_failure_reasons.csv`; the audit does not reinterpret them as solver failures.

## Representative curve limitation

The current aggregate artifact retains boundary/spinodal/CEP/crossover/grid CSVs and shard diagnostics, but not the raw rho-mu curve points used to draw S-shape figures. The nine regression anchors are therefore recorded as `not_available` in `tables/representative_curve_evidence.csv`; no historical curve is reused as a C2 plot.

## Author decision boundary

The next numerical action is not automatic. Before any rerun or policy change, review the regression geometry, the three retained CEP failures, and the failure-reason distribution. Do not relax tolerances, increase caps, promote phase-reference, or start RS production from this package alone.

## Files

- `manifest.json`: hashes, provenance and validation inventory
- `tables/`: derived gate, geometry, CEP, axis/reason, state/certificate, telemetry and evidence-gap tables
- `figures/`: diagnostic plots only; not formal publication figures
- `claim_ledger.csv`: supported, unsupported and unavailable claims
