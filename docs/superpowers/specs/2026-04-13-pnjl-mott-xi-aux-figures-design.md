# PNJL Mott xi auxiliary figures design

## Background and objective

The current analysis doc `docs/analysis/pnjl_mott/xi_dependence_analysis.md` already contains the core interpretation for xi dependence around the Mott transition, and `mode_ab` figures are available as merged comparative plots.

We now add an auxiliary figure set that improves evidence readability without re-running expensive scans. The user-approved execution order is:

1. Figure 1 (T_mott vs xi fit)
2. Figure 4 (order-parameter direct vs indirect effect)
3. Figure 2 (Gamma and Delta dual-axis)
4. Figure 5 (Taylor-like decomposition; only if feasible from existing data)

Figure 3 is not introduced as a new standalone figure. Instead, existing `mode_ab` outputs are reused with lightweight annotation copies to highlight crossing-point right shift and `Delta T_mott`.

## Scope

In scope:

1. Reuse existing CSV outputs and existing mode_ab PNGs.
2. Generate auxiliary figures under `docs/analysis/pnjl_mott/figures/` with agreed filenames.
3. Create annotated copies of existing mode_ab figures:
   - `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab/mott_mode_ab__M_K__xi3_annotated.png`
   - `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab/mott_mode_ab__M_pi__xi3_annotated.png`
4. Update `docs/analysis/pnjl_mott/xi_dependence_analysis.md` with an auxiliary-gallery section and references.
5. Normalize symbol convention in the document to:
   - `Delta = M_thr - M_mes`
   - plus an explicit note describing conversion from previous sign convention if it appears in historical text.

Out of scope:

- Re-running full scan workflows by default.
- Replacing existing plotting stack.
- Deleting or renaming historical outputs in `mode_a/` or `mode_b/`.

## Data sources and reuse policy

Primary inputs:

- `data/outputs/results/relaxtime/mott_phase/reference_100_300_fine/mott_phase_scan.csv`
- `data/outputs/results/relaxtime/mott_phase/reference_100_300_fine/mott_phase_derived.csv`
- `data/outputs/results/relaxtime/scan/gap_transport_scan_step5_muB0_xi-0p3to0p3.csv`
- Existing mode_ab figures in `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab/`

Policy:

1. Prefer direct reuse and post-processing from existing artifacts.
2. If Figure 5 cannot be derived from available columns, stop and report a minimal incremental computation plan before any new run.

## Deliverables and file contract

Auxiliary output directory:

- `docs/analysis/pnjl_mott/figures/`

Required figure files:

1. `fig1_tmott_vs_xi_fit.png`
2. `fig4_orderparam_direct_indirect.png`
3. `fig2_gamma_delta_dualaxis.png`
4. `fig5_taylor_decomposition.png` (conditional on data feasibility)

Mode_ab annotation copies:

1. `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab/mott_mode_ab__M_K__xi3_annotated.png`
2. `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine/mode_ab/mott_mode_ab__M_pi__xi3_annotated.png`

Doc update target:

- `docs/analysis/pnjl_mott/xi_dependence_analysis.md`

## Figure-by-figure design

### Figure 1: T_mott vs xi fit

Purpose:

- Quantify monotonic xi trend and curvature strength of `T_mott(xi)`.

Visual contract:

- Scatter points from existing `T_mott` samples.
- Linear and quadratic fit overlays.
- Legend includes fit equations and basic goodness indicators if available.

### Figure 4: direct vs indirect order-parameter effect

Purpose:

- Separate and compare:
  1) fixed-temperature snapshots (direct structural response),
  2) values sampled along `T_mott(xi)` trajectory (indirect transition-shift response).

Visual contract:

- Side-by-side panels or clear style separation in one canvas.
- Same y-units and explicit labels to prevent interpretation drift.

### Figure 2: Gamma and Delta dual-axis

Purpose:

- Jointly show broadening (`Gamma`) and threshold distance (`Delta`) around transition-sensitive region.

Visual contract:

- Left y-axis: `Gamma`.
- Right y-axis: `Delta = M_thr - M_mes`.
- Consistent xi color mapping with Figure 1 and Figure 4.

### Figure 3 handling (no standalone new figure)

Decision:

- No standalone Figure 3 generation.
- Reuse existing `mode_ab` K/pi figures and produce `_annotated` versions marking crossing points and right-shift relation.

### Figure 5: Taylor-like decomposition (conditional)

Purpose:

- Attribute xi dependence into interpretable decomposition terms.

Execution gate:

1. Attempt derivation from existing fields first.
2. If decomposition terms are not reconstructable, pause and present a minimal new-computation proposal (inputs, runtime estimate, expected outputs) for explicit approval.

## Architecture and implementation boundaries

1. Keep plotting logic in scripts under `scripts/analysis/` (new script or extension of existing analysis script), not under `tests/`.
2. Keep document rendering references relative to `docs/analysis/pnjl_mott/`.
3. Keep mode_ab annotation operation lightweight (image overlay/post-processing), independent from physics recomputation.

## Error handling and quality rules

1. Input CSV missing columns -> fail fast with explicit missing-column list.
2. Figure-level failure should report which figure failed and why; continue independent figures when safe.
3. If Figure 5 gate fails, mark as skipped with reason in both script log and markdown section.
4. All generated filenames must match the agreed contract exactly.

## Validation strategy

1. Deterministic local generation using existing artifacts.
2. Post-generation checks:
   - all expected files exist,
   - markdown references resolve to existing files,
   - symbol convention in text uses `Delta = M_thr - M_mes` consistently.
3. If tests are added, they should be lightweight and file-contract oriented.

## Acceptance criteria

1. Figures 1, 4, 2 are generated in agreed order with agreed names under `docs/analysis/pnjl_mott/figures/`.
2. Mode_ab annotated copies are generated with agreed `_annotated` names.
3. `xi_dependence_analysis.md` includes an auxiliary gallery section and references these assets.
4. No new scan run is performed unless Figure 5 proves impossible from existing data and explicit approval is obtained.
5. Sign convention is unified to `Delta = M_thr - M_mes`, with a clear conversion note for legacy wording.
