# Meson Density Literature Targets

This directory stores structured meson-density literature targets used by
validation-facing smoke checks.

Current role:

- preserve normalized target CSV assets under the formal validation data tree
- keep the matching `.meta.md` provenance beside each curve
- avoid treating these charged/freeze-out curves as regression truth

Current imported targets from `Blaschke:2019col`, Figure 4 right panel:

- `blaschke2019col_kminus_piminus_mu_pi_100_fig4_right.csv`
- `blaschke2019col_kminus_piminus_mu_pi_134p5_fig4_right.csv`
- `blaschke2019col_kplus_piplus_mu_pi_100_fig4_right_with_anomalous.csv`
- `blaschke2019col_kplus_piplus_mu_pi_100_fig4_right_no_anomalous.csv`

Supporting files:

- `manifest_wpd_import_report.csv`: normalized import report from the raw WPD export
- matching `.meta.md` files for each target

Interpretation note:

- these targets remain literature-result assets for workflow reproduction and
  data governance smoke checks
- they are not the first formal validation kernel target for PR98 follow-up
