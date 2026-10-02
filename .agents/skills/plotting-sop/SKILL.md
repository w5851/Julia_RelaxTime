---
name: plotting-sop
description: Apply the repository's current figure_production_v2 contract to new scientific plots, including frozen-input provenance, PNG review, and later vector delivery. Do not use for solver or numerical production changes.
---

# Project Plotting SOP

Use this skill when creating or reviewing a new paper-facing figure from
already frozen CSV/JSON results. The authoritative workflow is
`docs/guides/sop/workflows/figure_production.md`.

## Current contract

- New cases use `config/plotting/candidate_aps_v2.toml` for author review or
  `config/plotting/strict_aps_v2.toml` only when the input qualification and
  final-delivery scope justify strict status.
- The v1 profiles (`audit_v1`, `candidate_origin_like_v1`, and
  `strict_origin_like_v1`) are deprecated compatibility profiles. They may be
  loaded to reproduce historical artifacts, but never choose them for a new
  case.
- Reuse `scripts/plotting/plot_style.py`, `plot_manifest.py`,
  `plot_quality.py`, and `validate_plot_artifact.py`. A figure-family script
  owns data selection and physical labels; the shared layer owns the rendering
  and artifact contract.
- Do not call a solver, alter numerical CSV/JSON values, introduce hidden
  interpolation or smoothing, or modify the paper project as part of figure
  production.

## Two-stage delivery

1. Create a new sibling output directory. Refuse to overwrite an existing
   case. Record every input path, byte count, SHA-256, generator hash, Git
   commit, runtime, units, transform, series state, and mask/interpolation
   rule. For a retained historical case, validate code and plotting-contract
   hashes from the recorded full commit or source snapshot; validate current
   data, images, and manifests by their recorded bytes/hashes.
2. Generate the `png_review` stage first. It must contain PNG only and record
   `delivery_stage=png_review`, `manuscript_eligible=false`,
   `current_publication_layer=false`, and `vector_delivery_pending=true` in
   the manifest. Run all non-PDF quality checks: physical size, effective
   resolution, measured glyph height, line width, four-sided inward major and
   minor ticks, clipping, overlap, units, and output hashes.
3. Stop for author visual review. Review the PNG together with its manifest;
   do not promote it or update `publication_clean_current.json`.
4. After explicit acceptance, generate PDF from the same frozen inputs and
   the same plotting logic. Do not reselect data or change display semantics
   between stages. Run the PDF-only checks for one page, embedded non-Type-3
   fonts, no raster-wrapped curves, fixed physical dimensions, and matching
   hashes. Add EPS/PS only when the selected journal delivery route requires
   it.

## Visual requirements

- Use mathematical variable labels and put dimensional units in parentheses,
  such as `$\tau_u\;(\mathrm{fm})$` in the rendered label convention used by
  the figure family; never introduce square-bracket units in a new paper
  figure.
- Use inward ticks on top, bottom, left, and right. Linear axes use one
  minor tick between adjacent major ticks by default; log axes use the
  configured log minor locator. A denser linear minor-tick policy requires a
  case-level justification in the manifest.
- Draw at the intended physical insertion width. Measure actual glyph outline
  heights, not nominal em size; the v2 gate requires at least 2 mm for the
  measured capital/numeral glyphs. Dense PNG-only composites may use the
  explicitly documented review typography exception, with measured scripts
  at least 1.5 mm. This is not APS submission compliance and cannot promote
  the artifact or authorize vector delivery.
- Keep legends outside data axes unless a documented case layout and visual
  review justify a registered in-axes exception. Verify curve/landmark
  visibility, panel labels and clipping. Use color plus line style where
  curves must remain distinguishable in grayscale.
- Choose axis scales, ranges, legend hosts and precision for the scientific
  case. Tick labels must identify the actual values; display formatting does
  not establish numerical uncertainty.
- Keep figure-family layout details in its case documentation and renderer.
  The reviewed phase-guided v11 layout is documented in
  `docs/analysis/relaxtime/phase_guided_transport/plotting_case_contract.md`.
- Keep scientific phase labels and branch gaps tied to the input semantics;
  a plotting style change must not relabel a hadron/quark result or reconnect
  an audited discontinuity.

Historical manifests are evidence records, not requests to rerun the current
renderer. Use `validate_plot_artifact.py --snapshot --code-ref <full-sha>` for
the retained snapshot contract. A later change to the active SOP or generator
may require a new case or an explicit migration, but does not invalidate the
old image, input, or manifest bytes.

## Validation

Validate every new per-figure manifest and its actual output. Run public
contract/quality tests when shared plotting code, profiles or validators
change; run the affected figure-family tests when its logic changes.
Documentation and script-entrypoint checks follow changes to those contracts.
Generating a case with an existing renderer does not require every shared
test or Julia governance check.

The resulting package remains review-only even when all visual gates pass.
Formalization, manuscript eligibility, and current-layer replacement are
separate author-governed actions.
