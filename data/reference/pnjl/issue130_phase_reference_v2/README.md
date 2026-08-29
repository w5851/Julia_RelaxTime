# Issue #130 phase-reference v2

The public semantic layers are `strict`, `render` and `accepted`. The former v1 `derived` layer remains an internal build input and provenance record. `render` is the complete structured display package; `accepted` is the author-accepted default for downstream phase-map/analysis consumers. All layers are solver-free and `runtime_consumption=false`; the solver runtime continues to use strict certified rows with the explicit legacy fallback/rollback contract.
