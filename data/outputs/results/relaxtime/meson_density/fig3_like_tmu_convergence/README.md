# FIG3-like T-mu Convergence Audit

Scope: `phase_shift_gbu_reference`, `phase_display=fold_0_pi`, default equal-flavor chemical profile.

## Key Findings

- `nodes_40_vs_48_qmax12_omega10` max relative K/pi difference: `0.0643378`
- `nodes_32_vs_40_qmax12_omega10` max relative K/pi difference: `0.102861`
- `cutoff_8_vs_10_q32` max relative K/pi difference: `0.0611371`
- `cutoff_10_vs_q12omega10_q32` max relative K/pi difference: `0.120235`

## Production Choice

- Formal production uses `qmax=12`, `omega_max=10`, `q_nodes=48`, `omega_nodes=48`.
- Representative `40 -> 48` node sensitivity at this cutoff is below about 6.5% on the tested grid.
- Earlier `qmax=4`, `omega_max=4` is not adequate for FIG3-like production; cutoff differences versus 6 and above are large.

## Files

- `convergence_summary.csv`: pointwise comparisons.
- `q12w12/combined_meson_density_scan.csv`
- `q16w16/combined_meson_density_scan.csv`
- `q20w20/combined_meson_density_scan.csv`
- `q24w24/combined_meson_density_scan.csv`
- `q32w32/combined_meson_density_scan.csv`
- `q32w32_cutoff10/combined_meson_density_scan.csv`
- `q32w32_cutoff12_omega10/combined_meson_density_scan.csv`
- `q32w32_cutoff6/combined_meson_density_scan.csv`
- `q32w32_cutoff8/combined_meson_density_scan.csv`
- `q40w40_cutoff12_omega10/combined_meson_density_scan.csv`
- `q40w40_highmu_t160/combined_meson_density_scan.csv`
- `q48w48_cutoff12_omega10/combined_meson_density_scan.csv`
- `q8w8/combined_meson_density_scan.csv`
