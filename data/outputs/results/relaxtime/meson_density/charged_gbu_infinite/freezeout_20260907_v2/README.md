# Charged GBU smoke production result

## Adoption status

This directory is the retained ten-point charged infinite-thermal GBU result
used for smoke-production display and review. The current route configuration
authorizes `Models.run_charged_gbu_freezeout_scan` as the charged meson-density
smoke default with `production_default=true` and `default_tier="smoke"`.

The files in this directory are an immutable result package from the earlier
accepted run. Its original `manifest.json`, `run.json`, `ratios.csv`, density
tables, figures, source snapshot, and sidecar hashes are intentionally not
rewritten. The original run predates the default-marker adoption and therefore
retains `production_default=false` in its historical metadata. This README is
the adoption record; it is not a replacement manifest or a claim that the run
was re-executed after the marker change.

The current figure-side copies are at
`data/outputs/figures/relaxtime/meson_density/charged_gbu_infinite/freezeout_20260907_v2/`
with a standard `plot_manifest.json`. The original result-side PNG/PDF remain
inside this immutable bundle because they are covered by the original
146-entry output hash set. The source snapshot is retained in the external
`D:\\w\\jrt-ord` diagnostic archive and is not duplicated into the staged Git
result package.

## Physics scope

- background: `FixedMuBConservedCharges` quark-only BQS;
- constraints: `rho_Q/rho_B=0.4`, `rho_S=0`;
- meson feedback: disabled;
- route: infinite PNJL thermal difference with the two-line vacuum regulator,
  finite-q charged GBU shell integration, and explicit pion/kaon charged
  channels;
- coverage: 10 freezeout energies and 40 charged channel checks.

The observable is the GBU meson partial yield on the fixed quark-only
background. It is not a final-state hadron yield, a meson-feedback equilibrium
solution, or a formal numerical regression baseline.

## Provenance and checks

- original run status: `complete_research_curve_accepted`;
- original source head: `fa34a3975d858ca5c517f3490476178720f0d23e`;
- original manifest SHA-256:
  `2c180374b66ce5a0563a44be1f697341cdbb90910bc700e1be0ce6b926be2078`;
- the 146 output hashes recorded by the original manifest were checked against
  the retained files before adoption;
- the representative figure is `freezeout_ratios.png`; numeric ratios are in
  `ratios.csv` and per-energy channel evidence is retained in the same tree.

The default adoption does not promote this result to the repository's formal
baseline. Convergence, model-dependence, and experimental final-state mapping
remain separate review questions.
