# Branch Repair Evidence: 2026-06-20

This directory records the targeted repair for three branch-unstable
`FixedAsymmetricRho` points in
`trho_asymmetric_kplus_piplus_scan_v1`.

The repair used current `origin/main`
`36ef39d5cc57313d147a06175f9272260ec4c1f5`, where `trho_asymmetric`
equilibrium selection uses `pressure_max_all_attempts_multiseed`.

## Scope

Only these points were recomputed and copied into the parent production CSV:

| subdir | `T_MeV` | `rho_target` |
| --- | ---: | ---: |
| `T120_rho035` | `120` | `0.35` |
| `T130_rho080` | `130` | `0.80` |
| `T210_rho090` | `210` | `0.90` |

Each subdir contains a standalone `combined_meson_density_scan.csv` generated
with the same v1 production resolution settings:

- `p_num=8`, `t_num=4`, `max_iter=20`
- `stable_q_nodes=768`
- `q_nodes=192`, `omega_nodes=192`
- `qmax=4.0 fm^-1`, `omega_min=0.05 fm^-1`, `omega_max=4.0 fm^-1`
- `density_policy=x_min_cut`, `bose_x_min=1e-2`

This was not a full-grid rerun.
