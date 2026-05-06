# Meson Density Phase Point AD vs FD Workflow Probe

## Configuration

- scheme: `current`
- eta: `1.0e-6`
- fd_step: `1.0e-5`
- T_MeV: `208.0, 210.0, 212.0`
- xi: `0.0, 0.1`
- q: `0.0, 0.5, 1.0`
- omega: `0.3, 0.6, 1.0`
- mesons: `pi, K`

## Top Findings

- successful rows: `108`
- failed rows: `0`
- `dphase_formula_abs_diff` max = `7.815970093361102e-14` at `(T_MeV=208.0, xi=0.0, meson=pi, q=0.5, omega=1.0)`
- `|dphase_ad - dphase_fd|` max = `9.847695959092194` at `(T_MeV=208.0, xi=0.0, meson=K, q=0.0, omega=0.6)`
- `|dphase_ad - dphase_raw_fd|` max = `9.84758439322962` at `(T_MeV=208.0, xi=0.0, meson=K, q=0.0, omega=0.6)`
- `|dReD_domega - dReD_fd|` max = `191.57063132155196` at `(T_MeV=212.0, xi=0.1, meson=pi, q=0.5, omega=1.0)`
- `|dImD_domega - dImD_fd|` max = `361.69792112626374` at `(T_MeV=208.0, xi=0.0, meson=pi, q=0.5, omega=1.0)`

## Reading

1. `dphase_formula_abs_diff` should stay near machine precision if the AD path and the `Re/Im D` analytic derivative formula are internally consistent.
2. If `dphase_abs_diff` remains visibly nonzero while (1) is tiny, the main gap is not `atan` branch handling alone; it already exists at the `ReD/ImD` derivative level.
3. In the current implementation, that usually means AD is differentiating a frozen-topology integrand, while FD also samples how singularity intervals / segmentation move with `omega`.
