# Issue #130 RS `publication_clean_v5` figure-label review candidate

## Purpose and boundary

This package is a label-only display derivative of the author-accepted
`publication_clean_v4` figure layer.  It reuses the v4 point table, display
adjustments, phase gaps, and audit evidence byte-for-byte.  It changes no raw
CSV, production registry, solver output, or numerical display value.

The two requested rendering changes are:

1. The six relaxation-time y-axis labels now include the unit `[fm]`.
2. The first-order endpoint legend entries now read `chirally restored branch
   endpoint` and `chirally broken branch endpoint`.  The latter is wording for
   the quark transport branch and does not claim that hadronic transport was
   calculated.

The v4 package remains unchanged and remains the current formal publication
layer until this candidate is explicitly reviewed and accepted.  v5 is
`manuscript_eligible=false`, `solver_called=false`, and is not a convergence
or physical-model update.

## Provenance

- parent package: `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v4`
- parent package manifest SHA256: `e78f2c4b727b16a9294c79c896cd6733a7566df06fefb696166084e7357bc7eb`
- parent plot manifest SHA256: `90ae66d4478d5cfc699bfbf35dd7c9bffbea00da0268c8b13bb42068235779e6`
- source point table and audit tables: copied from v4 without content changes
- publication figures: 72
- raw/production data modified: false
- solver called for this derivative: false
- canonical numerical data modified: false

The exact old/new strings are recorded in
`tables/v5_display_label_map.csv`.  The v5 plot and package manifests record
the updated PNG hashes and the unchanged v4 parent hashes.

## Reproduction

    python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v5.py
    python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v5.py
    python -m pytest tests/unit/python/test_phase_guided_publication_clean_v5.py
