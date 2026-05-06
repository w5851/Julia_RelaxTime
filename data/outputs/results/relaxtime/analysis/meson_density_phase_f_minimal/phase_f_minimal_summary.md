# Phase F minimal summary

Input window: T = 208:2:220 MeV

Artifacts:

- density figure: `phase_f_minimal_n_densities.png`
- ratio figure: `phase_f_minimal_kpi_ratio.png`
- merged data: `phase_f_minimal_curves.csv`

Quick read:

1. `stable` remains the lowest-complexity baseline, with `K/pi` rising from about `0.57` to `0.64`.
2. `strict BW Stage2` is not close to `stable` in the current window; it enhances `n_K` strongly and keeps `K/pi` above `1`.
3. `current BU` stays much larger than `stable` in both `n_pi` and `n_K`, with `K/pi` around `1.37-1.40`.
4. `gbu reference` pulls the BU result back toward the `stable` scale and sits well below both `current BU` and `strict BW Stage2` in `K/pi`.

Endpoint values:

- 208 MeV:
  - stable: n_pi=0.34479, n_K=0.19562, K/pi=0.56737
  - strict BW: n_pi=0.20968, n_K=0.58336, K/pi=2.78209
  - current BU: n_pi=2.82408, n_K=3.86521, K/pi=1.36866
  - gbu reference: n_pi=0.26150, n_K=0.20272, K/pi=0.77520

- 220 MeV:
  - stable: n_pi=0.35847, n_K=0.22867, K/pi=0.63792
  - strict BW: n_pi=0.60417, n_K=1.03490, K/pi=1.71292
  - current BU: n_pi=2.79809, n_K=3.91433, K/pi=1.39893
  - gbu reference: n_pi=0.20833, n_K=0.18496, K/pi=0.88782
