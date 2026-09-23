# Issue #130 RS `publication_clean_v4` formal publication-clean layer

## Acceptance and scope

This is the author-accepted current publication-clean figure layer, promoted
from the reviewed v4 candidate on `2026-09-19T04:41:40.584741+00:00`.  It is a solver-free,
display-only derivative of the approved raw transport case.  The raw CSVs,
production registry, canonical data, and solver outputs are unchanged.

The v4 display derivative is `author_accepted_display_only` and is
`manuscript_eligible=true` for the disclosed display layer.  Source raw
numerical data remain `diagnostic_only` and `raw_manuscript_eligible=false`;
this decision does not certify a new solver result or convergence.  The local
high-rate gate was not run and is not represented as passed.

Manuscript eligibility covers the 72 v4 figures and the explicit raw/display
provenance in `tables/publication_clean_points.csv` and
`tables/v4_display_adjustment_map.csv`.  The mode-B composite mechanism
remains unassessed.  The mode-A `tau_sbar` endpoint remains a display
extrapolation, not a replacement branch solution or a basis for an exact jump
amplitude.  The first-order gap remains unfilled.

The eligibility decision is recorded in
`../publication_clean_v4_manuscript_eligibility_v1.json`.

## Accepted display rules

| scope | rule |
| --- | --- |
| mode-B composite curves | Linear interpolation between xi=0.35 and xi=0.37 at T=200 MeV, muB=0 MeV for zeta, sigma/T, and sigma. |
| mode-A tau_sbar endpoint | Log-linear extrapolation from the left-branch anchors xi=-0.02 and xi=-0.01 to xi=-0.003. |
| phase gate | Crossover-only interpolation; the first-order gap [-0.003, +0.003] remains a hard split. |

The endpoint raw value remains in `tables/publication_clean_points.csv` as
`raw_value`; the display value is retained separately as `clean_value`.  Do
not use the display endpoint for exact branch derivatives or jump amplitudes.

## Provenance

- promoted source: `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v4/`
- v4 source manifest SHA256: recorded in `publication_clean_formalization_v1.json`
- immutable intermediate parent: `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v3/`
- solver called for this layer: `false`
- raw/production data modified: `false`
- formal figure status: `author_accepted_formal_layout`
- numerical status: `diagnostic_only`

The v3 package is retained unchanged as the parent and intermediate evidence
snapshot.  The v4 package is the current publication layer; no lower-version
directory is overwritten or silently reinterpreted.

## Reproduction and audit

The v4 source build remains reproducible with:

    python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v4.py
    python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v4.py

The explicit author-acceptance record is
`docs/analysis/relaxtime/phase_guided_transport/publication_clean_formalization_v1.json`.
