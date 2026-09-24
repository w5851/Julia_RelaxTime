# Issue #130 RS `publication_clean_v5` formal publication-clean layer

## Acceptance and scope

This is the author-accepted current publication-clean figure layer, adopted
for the formal manuscript on `2026-09-24T10:38:42Z`.  It is a solver-free, label-only
display derivative of the accepted `publication_clean_v4` layer.  The raw CSVs,
production registry, canonical data, display values, and solver outputs are
unchanged.

The accepted rendering changes are:

1. The six relaxation-time y-axis labels include the unit `[fm]`.
2. The first-order endpoint legend entries read `chirally restored branch
   endpoint` and `chirally broken branch endpoint`.  The latter describes the
   displayed quark transport branch and does not claim that hadronic transport
   was calculated.

`manuscript_eligible=true` applies only to this disclosed v5 display layer.
The source raw numerical results remain `diagnostic_only` and
`raw_manuscript_eligible=false`.  The local high-rate convergence gate was not
run and is not represented as passed.

## Provenance

- parent display layer: `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v4`
- parent manuscript eligibility record: `docs/analysis/relaxtime/phase_guided_transport/publication_clean_v4_manuscript_eligibility_v1.json`
- publication figures: 72 (36 mode A and 36 mode B)
- raw/production data modified: false
- solver called for this derivative: false
- canonical numerical data modified: false
- label audit: `tables/v5_display_label_map.csv`
- raw/display provenance inherited from v4: `tables/publication_clean_points.csv` and `tables/v4_display_adjustment_map.csv`

The v4 parent remains retained unchanged as provenance.  No first-order gap is
filled, no branch solution is replaced, and no display label change is a
physical-model or convergence claim.

## Records and reproduction

The formalization record is
`../publication_clean_v5_formalization_v1.json`; the manuscript display
eligibility record is `../publication_clean_v5_manuscript_eligibility_v1.json`.

The candidate layer was generated with:

    python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v5.py
    python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v5.py
