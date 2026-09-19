# RS publication_clean_v4 formal publication figure layer

This directory records the byte-preserving publication figure mirror for
publication_clean_v4.

The v4 analysis package inherits v3, adjusts the three mode-B composite curves
at T=200, muB=0, xi=0.36, and applies one one-sided left-branch display
adjustment to mode-A tau_sbar at the audited first-order left endpoint.
The first-order gap remains split. Raw results, production registries, and
solver outputs are unchanged.

The layer is author-accepted as the current formal publication layout.
`manuscript_eligible=false` still applies to raw numerical claims: the
underlying layer remains solver-free and numerically diagnostic-only. The
mechanism review is retained as diagnostic context; it is not a propagator
regularization or a production convergence certificate.

Acceptance record:

    docs/analysis/relaxtime/phase_guided_transport/publication_clean_formalization_v1.json

Reproduction:

    python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v4.py
    python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v4.py
