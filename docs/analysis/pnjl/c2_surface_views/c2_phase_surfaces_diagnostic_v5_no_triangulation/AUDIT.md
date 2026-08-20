# C2 phase curves v5: no-triangulation audit

## Purpose

This package is a solver-free visualization audit of the fixed C2 reference
artifact. It tests whether the visible holes in the previous 3-D diagnostic
figure were introduced by `plot_trisurf` or are already present in the native
support. It does not create new numerical rows, reclassify an unresolved
point, or authorize phase-reference promotion.

## Input and provenance

- C2 source run: `31862752226`
- Calculation SHA: read from the immutable C2 reference manifest
- Source artifact: `stagec_density_v2_c2_20260813_run31862752226`
- Solver called: `false`
- Figure mode: `visual_lines`
- Triangulation: disabled globally

## Unresolved Maxwell semantics

`phase_grid_convergence.csv` is a convergence-diagnostic table, not the
Maxwell boundary table. An unresolved row can preserve refinement-level
metadata (`left`, `right`, `midpoint`) and/or geometry diagnostics
(`position_error_MeV`, `density_error`, `maxwell_area`). The bracket fields
are refinement-level coordinates, not a chemical-potential Maxwell crossing
bracket. These fields describe the last evaluated comparison, but they do not
make the row a valid Maxwell boundary point.

Only a finite, converged row in `boundary_*.csv` is treated as native Maxwell
support. Therefore an unresolved grid record may have useful residuals while
having no plotted Maxwell point. The generated
`tables/grid_unresolved_diagnostics.csv` records this distinction for every
unresolved row, including the bracket semantics and whether a matching
converged boundary row exists.

For the fixed C2 artifact, the unresolved population is 5,424 rows:

- 5,337 have all six candidate-diagnostic fields present;
- 87 have only a partial set of those fields;
- the rho/temperature rows cover 3,485 unique `(xi,T)` coordinates; 612 of
  those coordinates have no converged boundary row;
- 4,253 share their `(xi,T)` coordinate with some converged boundary row
  (typically a different refinement level), while 1,171 have no matching
  converged boundary row;
- none of the unresolved rows are silently promoted to a boundary point.

`boundary_row_present` in the audit table is therefore coordinate-level
provenance, not a claim that the unresolved level itself passed. The native
boundary file contains 6,886 rows and every one is `converged=true`; an
unresolved level is never added to that file merely because it has a finite
area or a finite crossing bracket.

For example, the `xi=0.5`, `T=17 MeV` rho records retain finite level-3/4
geometry diagnostics but are marked `rho_geometry_not_converged` and
`hybrid_stage_c_not_converged`; the absence of a matching boundary row is
intentional under the production contract.

## Rendering contract

- Maxwell and crossover are drawn only from native, ordered support rows.
- Adjacent points are connected only when their native coordinate gap is
  within the explicit plotting limit; larger gaps remain visibly open.
- No triangle mesh, surface fill, interpolation, synthetic point, or gap
  connector is generated.
- CEP information remains a bracket projection and midpoint marker; the
  midpoint is not a certified CEP.
- Response peaks above the per-`xi` CEP chemical-potential proxy remain in
  the audit tables but are not drawn as physical crossover.

The resulting figure therefore separates two cases: a hole that disappears
when triangulation is removed was a display artifact, while a hole that
remains in v5 reflects missing native support or physical filtering.

## Boundary

This package is diagnostic-only. It must not be used to promote unresolved
geometry/interpolation records, infer a continuous phase surface across a
gap, or replace the C2 convergence and author-review gates.
