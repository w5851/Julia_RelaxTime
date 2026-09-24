# RS publication_clean_v5 figure-layer review candidate

This directory records the byte-preserving public mirror for the v5 label-only
derivative of publication_clean_v4.  It adds `[fm]` to relaxation-time y-axis
labels and replaces `quark`/`hadron` endpoint wording with chirally restored /
chirally broken branch endpoint wording.  Numerical values, raw data, audit
tables, phase gaps, and solver outputs are unchanged.

The layer is for author review only and remains `manuscript_eligible=false`.

Reproduction:

    python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v5.py
    python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v5.py
