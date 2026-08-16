# first_canonical_v1_p128_xi005_validated_anchored_prod_v1 Convergence Evidence

Shared convergence and artifact audit records for the promoted p128 phase-guided transport production case.

## Verdict

`production-grade`

The promoted result is based on the `128/20/36/24/560` parameter set:

- `tau_p_nodes = 128`
- `tau_angle_nodes = 20`
- `tau_phi_nodes = 36`
- `tau_n_sigma_points = 24`
- `sigma_grid_n = 560`

The adjacent high-precision reference is `104/18/30/20/440`.

## Files

- `p104_vs_p128_convergence_gate.csv`: pointwise convergence comparison
- `p104_vs_p128_convergence_gate.summary.json`: machine-readable convergence summary
- `p104_vs_p128_convergence_gate.md`: human-readable convergence summary
- `audit_summary_refined.json`: merged artifact audit summary
- `download_manifest.csv`: p128 shard labels and GitHub Actions run ids
- `merged_manifest.json`: merged mode-level row counts and source directories

## Gate Summary

- matched points: `378`
- missing p104 points: `0`
- missing p128 points: `0`
- overall max relative difference: `0.009426486578504083`
- core max relative difference: `0.006176952435828459`
- anti-particle max relative difference: `0.009426486578504083`
- overall count above 1%: `0`
- overall count above 0.5%: `12`

Worst point:

- `mode_a`, `muB=450 MeV`, `alpha_T=1.1`, `xi=0.35`, `tau_ubar/tau_dbar`, relative difference `0.9426486578504083%`

## p128 Shards

The p128 production artifact was generated as 18 GitHub Actions shards at commit `c3f25dde5c9e62e7a3ebcd8ed50aef2b53377c0f`.

| shard | run id |
| --- | --- |
| `a_muB0_alpha1p0` | `27246554475` |
| `a_muB0_alpha1p1` | `27246561257` |
| `a_muB0_alpha1p2` | `27246567405` |
| `a_muB450_alpha1p0` | `27246574295` |
| `a_muB450_alpha1p1` | `27246582982` |
| `a_muB450_alpha1p2` | `27246592047` |
| `a_muB900_alpha1p0` | `27246598301` |
| `a_muB900_alpha1p1` | `27246605064` |
| `a_muB900_alpha1p2` | `27246611493` |
| `b_T120_muB0` | `27246617302` |
| `b_T120_muB450` | `27246622499` |
| `b_T120_muB900` | `27246627637` |
| `b_T160_muB0` | `27246634399` |
| `b_T160_muB450` | `27246639390` |
| `b_T160_muB900` | `27246646432` |
| `b_T200_muB0` | `27246651646` |
| `b_T200_muB450` | `27246657527` |
| `b_T200_muB900` | `27246669072` |

The p104 comparison artifacts were also generated as 18 shards. Their run ids are preserved in the local audit record and summarized in `p104_vs_p128_convergence_gate.*`.
