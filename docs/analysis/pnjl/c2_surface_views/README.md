# C2 Surface Views

本系列的统一只读 namespace 索引为 [`../phase_surface_series/`](../phase_surface_series/)。本目录仍是历史 source package；其 manifest、checksum、execution log、生成时 provenance 和原始路径均保持不变。namespace snapshot 的逐文件 SHA-256 对照见 [`../phase_surface_series/series_manifest.json`](../phase_surface_series/series_manifest.json)。

This group contains the versioned diagnostic surface line for C2:

- `c2_phase_surfaces_diagnostic_v1/` and `c2_phase_surfaces_diagnostic_v2/`: chemical-potential convention corrections.
- `c2_phase_surfaces_diagnostic_v3/`: crossover physical filtering.
- `c2_phase_surfaces_diagnostic_v4_visual_closed/`: visual-closure rendering.
- `c2_phase_surfaces_diagnostic_v4_visual_closed_display16/`: the v4 display16 variant with its own audit and manifest.
- `c2_phase_surfaces_diagnostic_v5_no_triangulation/`: native-support/no-triangulation rendering.
- `c2_phase_surfaces_diagnostic_v6_crossover_overlay/`: v5 postprocessed baseline with the solver-free-replayed crossover endpoint overlay; Maxwell, spinodal and CEP evidence are inherited unchanged.
- `c2_phase_surfaces_diagnostic_v7_crossover_derived/`: v6-based solver-free derived crossover layer. Same-`xi` gaps and adjacent-`xi` common-support gaps are piecewise-linearly filled with explicit `interpolated_noncertified` provenance; CEP boundary rows use an `estimated_midpoint` boundary estimate. Maxwell remains native v6 support only.

The eight cases share the C2/v5 baseline source run `31862752226`, but each preserves its own semantic verdict, tables, decision record, and manifest. v6 additionally records the endpoint expansion numerical/replay provenance; v7 records the derived-layer source keys and coverage mask. They are `diagnostic_only` and do not promote unresolved geometry, CEP midpoints, endpoint overlays, or interpolated values to phase-reference truth.

The move changes only the repository namespace and current references. Per-case manifests retain generation-time external provenance paths, and frozen audit packages retain their historical paths and hashes. A pre-move check found existing `manifest.json` `output_files` hash mismatches in all six cases; this namespace-only batch does not rewrite or normalize those evidence records. Treat that as a separate provenance-audit follow-up.
