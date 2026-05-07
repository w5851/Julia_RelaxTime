# Plot Review Case: freezeout_kminus_piminus_mu_pi_100

Current role:

- plot-review asset for manual trend inspection
- not an external validation gate
- not a regression truth source by itself

Inputs:

- workflow source: `data/outputs/results/relaxtime/meson_density/freezeout_validation/blaschke2019col_kminus_piminus_mu_pi_100_phase_shift_gbu_default/workflow_scan.csv`
- comparison source: `data/outputs/results/relaxtime/meson_density/freezeout_validation/blaschke2019col_kminus_piminus_mu_pi_100_phase_shift_gbu_default/comparison_vs_target.csv`

Outputs:

- `workflow_scan.csv`
- `comparison_vs_target.csv`
- `plot_review_comparison.csv`
- `overlay_kminus_piminus_mu_pi_100.png`
- `residual_kminus_piminus_mu_pi_100.png`

Manual review checklist:

- compare target/model monotonic rise
- compare whether the workflow remains systematically below the literature target
- compare whether no nonphysical reversals or spikes appear
- compare whether the relative-difference band shape changes after code updates

Summary metrics:

- points: `48`
- max abs diff: `0.147787`
- max rel diff (finite only): `0.873012`

Interpretation note:

- This directory supports regression fallback and manual plot review only.
- It does not promote the literature curve to `tests/validation/`.
