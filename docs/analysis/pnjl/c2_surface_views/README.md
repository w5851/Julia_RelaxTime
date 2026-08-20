# C2 Surface Views

This group contains the versioned diagnostic surface line for C2:

- `c2_phase_surfaces_diagnostic_v1/` and `c2_phase_surfaces_diagnostic_v2/`: chemical-potential convention corrections.
- `c2_phase_surfaces_diagnostic_v3/`: crossover physical filtering.
- `c2_phase_surfaces_diagnostic_v4_visual_closed/`: visual-closure rendering.
- `c2_phase_surfaces_diagnostic_v4_visual_closed_display16/`: the v4 display16 variant with its own audit and manifest.
- `c2_phase_surfaces_diagnostic_v5_no_triangulation/`: native-support/no-triangulation rendering.

The six cases share C2 source run `31862752226`, but each preserves its own semantic verdict, tables, decision record, and manifest. They are `diagnostic_only` and do not promote unresolved geometry, CEP midpoints, or visual closure to phase-reference truth.

The move changes only the repository namespace and current references. Per-case manifests retain generation-time external provenance paths, and frozen audit packages retain their historical paths and hashes. A pre-move check found existing `manifest.json` `output_files` hash mismatches in all six cases; this namespace-only batch does not rewrite or normalize those evidence records. Treat that as a separate provenance-audit follow-up.
