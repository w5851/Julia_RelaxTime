# Issue #130 phase-reference layers v1

This package freezes a strict input layer, a uniform-xi derived layer and a structured render layer. The strict layer contains immutable C2/v6 support plus the 276-target Maxwell expansion; the derived layer performs only local common-support interpolation; the render layer is a no-triangulation visual projection. All layers keep `reference_write=false` and remain diagnostic/author-review candidates.
