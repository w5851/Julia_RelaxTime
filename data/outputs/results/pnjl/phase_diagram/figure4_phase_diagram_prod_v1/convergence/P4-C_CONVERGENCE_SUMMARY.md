# P4-C PNJL Phase Reference Convergence Summary

Date: 2026-07-05

## Verdict

`production-grade`

This verdict applies to the PNJL phase reference data generated with the C1 full-grid numerical policy:

- `p_num=24`
- `t_num=8`
- `iterations=80`
- `T_min=60`, `T_max=240`, `T_step=5`
- `rho_min=0.0`, `rho_max=4.0`, `rho_step=0.05`
- `xi=-0.5:0.05:0.5`

The verdict is supported by the C2 anchor comparison at `xi=-0.5,-0.25,0.0,0.25,0.5` using `p_num=32`, `t_num=12`, and `iterations=100`.

## GitHub Actions Runs

- C1 full-grid candidate: https://github.com/w5851/Julia_RelaxTime/actions/runs/28736257287
- C2 refined anchors: https://github.com/w5851/Julia_RelaxTime/actions/runs/28741709519
- head SHA: `aca739f6bac4565dddef888abf9588b8c0fa583f`

## Downloaded Artifacts

- C1 artifact: `action_artifacts/c1_p24t8/`
- C2 artifact: `action_artifacts/c2_p32t12_anchors/`
- C1 vs C2 comparison: `c1_vs_c2_anchor_comparison/`

## C1 Inventory

- boundary rows: 304
- CEP rows: 21
- spinodal rows: 304
- crossover rows: 336
- crossover converged rows: 336 / 336
- xi coverage: `-0.5:0.05:0.5` (21 values)
- NaN/Inf count: 0 across boundary, CEP, spinodal, and crossover CSV files
- negative physical values outside `xi`: none detected in the downloaded CSV audit

## C2 Inventory

- boundary rows: 72
- CEP rows: 5
- spinodal rows: 72
- crossover rows: 80
- crossover converged rows: 80 / 80
- xi coverage: `[-0.5,-0.25,0.0,0.25,0.5]`
- NaN/Inf count: 0 across boundary, CEP, spinodal, and crossover CSV files
- negative physical values outside `xi`: none detected in the downloaded CSV audit

## C1 vs C2 Anchor Comparison

Machine-readable outputs:

- `c1_vs_c2_anchor_comparison/phase_reference_artifact_inventory.csv`
- `c1_vs_c2_anchor_comparison/phase_reference_convergence_comparison.csv`
- `c1_vs_c2_anchor_comparison/phase_reference_convergence_summary.json`
- `c1_vs_c2_anchor_comparison/phase_reference_convergence_summary.md`

Key maximum differences:

| artifact | metric | max abs diff | max rel diff |
|---|---:|---:|---:|
| boundary | `mu_transition_MeV` | 0.2613591590 | 7.180345729e-4 |
| boundary | `rho_hadron` | 0.0057633117 | 1.998829794e-2 |
| boundary | `rho_quark` | 0.0399753263 | 1.601860208e-2 |
| CEP | `T_CEP_MeV` | 0.0292968750 | 2.737850787e-4 |
| CEP | `muB_CEP_MeV` | 0.0759663347 | 7.269334828e-5 |
| crossover | `T_crossover_MeV` | 0.0304144661 | 2.842399250e-4 |
| crossover | `mu_MeV` | 0.0253221116 | 7.269334828e-5 |
| spinodals | `mu_spinodal_hadron_MeV` | 0.0716456155 | 1.903876033e-4 |
| spinodals | `mu_spinodal_quark_MeV` | 1.4341450249 | 4.188879452e-3 |
| spinodals | `rho_spinodal_hadron` | 0.0058768620 | 1.085970101e-2 |
| spinodals | `rho_spinodal_quark` | 0.0492897256 | 2.781935082e-2 |

## Residual Risks

- C2 is an anchor-grid refined check, not a full `xi=-0.5:0.05:0.5` refined rerun.
- Crossover `rho` and derivative diagnostics show larger relative differences at some high-order points; these are retained for audit but are not the primary `T-mu_B` phase-line observables.
- Spinodal data pass the anchor convergence check but should remain visually optional in the default figure because it can add density and branch-interpretation complexity.

## Decision

C1 is accepted as the production-grade full-grid phase reference source for P4-D, with the residual risks above recorded for downstream audit and figure design.
