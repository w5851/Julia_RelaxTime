# RS publication_clean_v4 formal publication figure layer

This directory records the byte-preserving publication figure mirror for
publication_clean_v4.

The v4 analysis package inherits v3, adjusts the three mode-B composite curves
at T=200, muB=0, xi=0.36, and applies one one-sided left-branch display
adjustment to mode-A tau_sbar at the audited first-order left endpoint.
The first-order gap remains split. Raw results, production registries, and
solver outputs are unchanged.

The layer is author-accepted and `manuscript_eligible=true` for display
use with raw/display values and local adjustments disclosed.  Source raw
numerical data remain `diagnostic_only` and are not promoted by this decision.
No local high-rate convergence gate was run; the figure status does not imply
that such a gate passed.

The mode-B composite adjustment has no separate mechanism verdict in the v4
record.  The `tau_sbar` endpoint remains a display extrapolation and cannot be
used as a replacement branch solution or an exact jump amplitude.

Eligibility record:

    docs/analysis/relaxtime/phase_guided_transport/publication_clean_v4_manuscript_eligibility_v1.json

Acceptance record:

    docs/analysis/relaxtime/phase_guided_transport/publication_clean_formalization_v1.json

Reproduction:

    python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v4.py
    python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v4.py
