# publication_clean_v10 PNG review layer

v10 is a solver-free display derivative of the accepted v5 point table. It
preserves all v5 CSV values, display-adjustment records, phase gaps, raw
provenance, and the current publication pointer. It emits 72 single PNG
charts (36 mode A and 36 mode B) plus two native composite PNG charts.

The layout follows a thesis-like case exception to the default external
legend policy: each single chart has one in-axes legend. Each composite has
the alpha_T key in (a) and the restored/broken marker key under a shared
1st-order heading in (c), once each. Panel labels are consistently above the
upper-left frame; no labels are moved onto low-valued curves.
Linear panels use three labeled x-axis major ticks, one linear minor tick
between adjacent major ticks, and a MaxNLocator linear y-axis. Figure 1's
complete muB=900 MeV column uses log-y; other positive panels with
data_max/data_min >= 20 use local log-y. The v10 manifest records the parent
scale and the v10 display scale for each panel. In-panel legends are checked
against rendered curves and endpoint markers, with zero intersections. The
sigma/T major ticks are exact at three-decimal precision with leading zeros.
Composite labels/titles/ticks are 11 pt; legend text is 10.5 pt. Main text
glyphs exceed 2 mm, but math scripts near 1.69 mm require the documented
PNG-review-only typography exception, not submission authorization.

This is a PNG-only author-review stage. It remains
`manuscript_eligible=false`, `current_publication_layer=false`, and
`vector_delivery_pending=true`; it does not replace publication_clean_v5.
No solver, high-rate gate, numerical production, paper-project edit, or new
display value was used.

Parent v5 manifest SHA256: `6f31a534bff0f75701091dbf1565614eb2e17290faa6cc2e0451c4703dc8fcc8`
Style manifest: `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v10_png_review/v10_display_style_manifest.json`
Caption handoff: `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v10_png_review/caption_handoff.md`

Reproduce in an absent sibling directory:

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v10.py --png-review
```
