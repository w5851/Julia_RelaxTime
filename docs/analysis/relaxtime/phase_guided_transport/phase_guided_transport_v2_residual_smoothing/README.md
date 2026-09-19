# `publication_clean_v2` residual-smoothing recipe

This directory contains the explicit, display-only recipe consumed by
`build_phase_guided_publication_clean_v2_residual_smoothed.py`.

The recipe is not a numerical result and does not authorize a solver rerun. Each
row names a curve-local anchor interval, the observables to display with a
shape-preserving local trend, and the blend strength. First-order endpoint gaps
remain hard barriers; a recipe row cannot fill or cross them.

Interpretation boundaries:

- `journal/main` rows are candidates for author review of the main figure layer,
  not manuscript eligibility or production promotion.
- `supplemental_author_review` is intentionally outside the main-figure claim
  set until the author accepts the visual change.
- The `mode_b / T=200 MeV / muB=900 MeV / xi=-0.10` row includes a
  `tau_sbar` channel-rate candidate whose local high-rate convergence gate is
  still empty. Smoothing that row changes only its display value; it does not
  resolve the mechanism or certify convergence.

The generated artifact preserves `raw_value` and the prior v2 `clean_value` and
records the final display value and point-level residual audit separately.
