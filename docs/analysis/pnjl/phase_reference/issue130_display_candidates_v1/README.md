# Issue #130 display-candidate provenance

This evidence group stores the non-figure provenance for the Figure 4 v3, v4,
and v5 display-only candidates. The corresponding PNG/SVG files and
`plot_manifest.json` remain under `data/outputs/figures/`; this separation is
required by the repository data-output path guard.

The candidates are style-only renders. They do not modify strict/derived CSV
data, call the solver, or switch the runtime reference consumer.
