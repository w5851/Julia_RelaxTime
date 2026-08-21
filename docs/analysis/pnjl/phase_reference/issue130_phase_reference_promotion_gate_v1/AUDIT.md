# Issue #130 phase-reference promotion gate v1

- Verdict: `promotion_candidate`
- Calculation SHA: `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- Numerical source run: `32354095831`; aggregate replay: `32451053476`
- Strict Maxwell unresolved rows retained: `4071`
- This gate is solver-free. It does not write `data/reference` and does not switch runtime consumers.
- A passing gate authorizes preparation of a separate versioned import PR only.

## Boundaries

- `strict_reference_v1` remains the evidence/source layer; unresolved geometry status is preserved.
- `derived_reference_v1` is an explicit downstream derived layer; interpolated rows are not certificates.
- `phase_surface_render_v1` is a reproducible render/data projection and cannot change physical status.
