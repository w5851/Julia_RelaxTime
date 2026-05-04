# Meson Density Literature Targets

This directory stores manually confirmed literature targets for meson-density / `K/pi` validation.

Suggested naming:

- `paperid_target_label.csv`
- `paperid_target_label.meta.md`

Recommended CSV schema:

```csv
x_value,y_value
```

Recommended meta fields:

1. source paper
2. figure number
3. panel
4. curve label
5. physical quantity
6. scenario qualifiers (`charged`, `freeze-out`, `mu_pi`, `mu_s`, `generalized_BU`, etc.)

Current imported targets from `Blaschke:2019col`, Figure 4 right panel:

- `blaschke2019col_kminus_piminus_mu_pi_100_fig4_right.csv`
- `blaschke2019col_kminus_piminus_mu_pi_134p5_fig4_right.csv`
- `blaschke2019col_kplus_piplus_mu_pi_100_fig4_right_with_anomalous.csv`
- `blaschke2019col_kplus_piplus_mu_pi_100_fig4_right_no_anomalous.csv`

Supporting files:

- `manifest_wpd_import_report.csv`: normalization/import report from the raw WPD export
- matching `.meta.md` files for each target

Current recommended first minimal validation target:

- `Blaschke:2019col`
- Figure 4 right panel
- `K^- / pi^-`
- `mu_pi = 100 MeV`

Important interpretation note:

- The `K^+ / pi^+` result in this figure is not a single line; it should be treated as a band-like result layer composed of:
  - `with_anomalous`
  - `no_anomalous`
- The right-panel scan path is currently understood only qualitatively as a stitched phenomenological path involving phase-transition and constant-temperature segments; it is not yet reconstructed as an exact parametric workflow target.
