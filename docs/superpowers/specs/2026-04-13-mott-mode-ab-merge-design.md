# Mott phase mode_ab merged comparison design

## Background and goal

Current plotting flow outputs:

- `mode_a`: one figure per `xi`, each figure contains all 6 observables.
- `mode_b`: one figure per observable, each figure contains all `xi` curves.

For visual comparison, we need a merged view that reduces clutter from the previous 5x6 combination. The new target is two figures total, each focused on one meson family and only three `xi` values.

## Scope

In scope:

1. Add a new output mode folder `mode_ab` in the existing Julia entry script.
2. Generate exactly two figures:
   - `mott_mode_ab__M_K__xi3.png`
   - `mott_mode_ab__M_pi__xi3.png`
3. Keep only `xi = -0.3, 0.0, 0.3`.
4. Use these observable triplets:
   - K figure: `M_K`, `Gamma_K`, `M_u_plus_M_s`
   - pi figure: `M_pi`, `Gamma_pi`, `M_u_plus_M_d`
5. Reuse existing derived CSV; no scan/data regeneration.
6. Preserve existing `mode_a` and `mode_b` behavior.

Out of scope:

- Changing scan configuration or recomputing physics data.
- Replacing the Python plotting backend.
- Renaming or removing existing `mode_a`/`mode_b` outputs.

## Architecture and data flow

1. Julia entry script `scripts/relaxtime/run_mott_phase_plot_modes.jl` remains orchestration layer.
2. Python script `scripts/plot_scan_csv.py` remains rendering layer.
3. New `_mode_ab(...)` function in Julia script will:
   - call Python plotting twice (for K group and pi group), using filtered `xi` via `--group xi` behavior consistent with existing mode flow,
   - stage outputs in temp directory,
   - copy/rename final figures into `mode_ab/`.
4. Existing `_mode_a(...)` and `_mode_b(...)` remain unchanged.

## Detailed plotting contract

### Input

- CSV: `data/outputs/results/relaxtime/mott_phase/reference_100_300_fine/mott_phase_derived.csv`

### Xi filtering

- Include only `-0.3`, `0.0`, `0.3`.
- Exclude `-0.15` and `0.15`.

### Figure grouping

- K figure uses y-list: `M_K,Gamma_K,M_u_plus_M_s`
- pi figure uses y-list: `M_pi,Gamma_pi,M_u_plus_M_d`

### Output path

- Root: `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine`
- New folder: `mode_ab`
- Files:
  - `mode_ab/mott_mode_ab__M_K__xi3.png`
  - `mode_ab/mott_mode_ab__M_pi__xi3.png`

## Error handling and fallback behavior

1. If required source PNG is not found in expected temp path, fallback to first PNG in directory (same strategy as existing `_mode_a`).
2. If no PNG exists after plotting, skip copy for that figure and keep script behavior non-destructive for other modes.
3. Input and script path validation remain unchanged (`ArgumentError` with explicit path).

## Test strategy

Update integration smoke test `tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl`:

1. Keep existing assertions for `mode_a` and `mode_b`.
2. Add assertions for:
   - `mode_ab/mott_mode_ab__M_K__xi3.png`
   - `mode_ab/mott_mode_ab__M_pi__xi3.png`

This verifies merged output generation without broadening test scope.

## Risks and mitigations

1. Risk: filename pattern differences from Python output internals.
   - Mitigation: retain fallback-first-PNG copy behavior and explicit destination names.
2. Risk: float matching for `xi=0.0`.
   - Mitigation: use tolerant comparison (`isapprox`) against target `xi` list.
3. Risk: regressions to existing modes.
   - Mitigation: preserve existing code paths and keep smoke assertions for old outputs.

## Acceptance criteria

1. Running existing plot mode script against current derived CSV creates old outputs plus new `mode_ab` outputs.
2. Only two new merged figures are added in `mode_ab` with agreed names.
3. Merged figures only include `xi=-0.3,0,0.3`.
4. K figure contains `M_K/Gamma_K/M_u_plus_M_s`; pi figure contains `M_pi/Gamma_pi/M_u_plus_M_d`.
5. Integration smoke test passes with new assertions.
