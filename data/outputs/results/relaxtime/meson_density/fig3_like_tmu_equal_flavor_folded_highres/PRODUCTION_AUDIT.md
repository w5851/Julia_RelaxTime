# FIG3-like Combined Meson Density Production Audit

## Scope

- Path strategy: `tmu`, multiple fixed `mu_q` temperature scans.
- Chemical profile: `default`, so `mu_u=mu_d=mu_s=mu_q` and `mu_pi=mu_K=0`.
- Phase policy: `real_axis_mode=pv_b0_eta0`, `phase_convention=arg_inverse_propagator`, `phase_display=fold_0_pi`.
- Numerical settings: `qmax=12`, `omega_max=10`, `q_nodes=48`, `omega_nodes=48`, `stable_q_nodes=256`.
- Grid: `T=120:10:240 MeV`, `mu_q=0:50:300 MeV`.

## Outputs

- `combined_meson_density_scan.csv`
- `combined_meson_density_scan.svg`
- `combined_meson_density_fig3_like.png`
- `run.stdout.log`, `run.stderr.log`

## Validation Summary

- Rows: `364`
- Status `ok`: `364`
- Negative density rows: `0`

## Phase-Display Diagnostic

- Representative point: `T=200 MeV`, `mu_q=0`, `phase_shift_gbu_reference`, `real_axis_mode=pv_b0_eta0`, `phase_convention=arg_inverse_propagator`.
- Legacy low-cutoff window `qmax=4`, `omega_max=4`, `16x16`: unwrapped phase gives `n_pi=-0.04425289187141765`, phase range `[-4.3886973001435745, 1.2216897319175892]`, negative phase samples `110/256`.
- Applying `phase_display=fold_0_pi` on the same samples gives `n_pi=0.11867743440016269`, folded phase range `[0.0, pi]`, out-of-range samples `0/256`.
- Production cutoff check `qmax=12`, `omega_max=10`, `16x16`: folded phase range remains `[0.0, pi]`, out-of-range samples `0/256`.
- Interpretation: negative `n_pi` is caused by feeding the unwrapped negative phase branch into the GBU weight; it is not an anomalous-mode effect in the equal-flavor chemical-potential path.

## K/pi Range By Regime

| regime | min | max | mean |
| --- | ---: | ---: | ---: |
| `phase_shift_current` | 1.31201 | 1.50807 | 1.38691 |
| `phase_shift_gbu_reference` | 0.689507 | 1.39083 | 0.891437 |
| `stable` | 0.196984 | 1.13997 | 0.525072 |
| `strict_bw_stage1` | 0.117696 | 1.16172 | 0.518037 |

## Convergence Boundary

- Convergence evidence is stored in `../fig3_like_tmu_convergence/`.
- Representative `q_nodes/omega_nodes 40 -> 48` sensitivity at `qmax=12`, `omega_max=10` is below about 6.5% on the tested grid.
- Earlier `qmax=4`, `omega_max=4` is not adequate for this FIG3-like production.
- This remains a mainline diagnostic production, not a claim that BU2020 FIG3 density-layer prescriptions are uniquely specified by the paper.
