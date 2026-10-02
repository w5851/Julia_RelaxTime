# publication_clean_v11 PNG review layer

v11 reuses frozen v5 data and v10's reviewed physical layout. Figure 1 is the
12-panel relaxation-time composite, now uniformly log-y with ordinary numeric
ticks. Figure 2 is the nine-panel transport composite, still linear; eta/s and
zeta/s share two decimal places within each row and sigma/T shares three.
Matching mode-A tau singles use log-y. The First-order endpoint key in (c)
is left-aligned and the same 10.5 pt size as the parameter key in (a).
Caption scope separates the common parameter key from muB900.0 alpha1.0
endpoint-bearing curves. Each panel retains its independent y range.

The package contains 72 single and two composite 600 dpi PNGs. It remains
manuscript_eligible=false, current_publication_layer=false, and
vector_delivery_pending=true. The v10 dense-composite typography review
exception remains explicit; this is not full APS submission compliance.
No v5 values, phase gaps, raw provenance, current pointer, prior v10 file,
solver, numerical gate, smoothing, PDF, or paper-project file was changed.

Manifest and caption: `manifest.json`, `v11_display_style_manifest.json`,
`caption_handoff.md`. Retained v10 artifact hashes are included as inputs.
This independent figure task does not replace the primary task-ledger track.

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v11.py --png-review
```
