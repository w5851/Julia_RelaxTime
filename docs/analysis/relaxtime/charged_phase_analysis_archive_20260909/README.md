# Charged phase analysis archive index

This directory is a small provenance index for the charged RPA/GBU freezeout
diagnostics retained during the PR #311 review. It is not a copy of the bulk
numeric archive and it does not promote any diagnostic result to a formal
baseline.

## Retained layers

- The accepted ten-point research result remains in
  `data/outputs/results/relaxtime/meson_density/charged_gbu_infinite/freezeout_20260907_v2/`.
  The route configuration is now the charged infinite-thermal GBU smoke
  default (`production_default=true`, `default_tier=smoke`). The immutable run
  manifest still says `production_default=false`; that is the identity of the
  earlier run, not the current entrypoint configuration.
- The finite-q/q=0 comparison and historical freezeout overlay remain in the
  analysis layer at
  `data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/`.
  Their own manifests all set `production_authorized=false`; the historical
  overlay is explicitly a context-only, unmatched projection.
- The analysis files are intentionally still untracked bulk evidence. This
  commit records their paths and hashes, but does not stage the PNG/CSV trees
  themselves. The source snapshots embedded in two analysis directories are
  also recorded only by count, size, and aggregate hash.

## Scope of the light index

The detailed inventory covers 32 evidence files (1,796,532 bytes):

- `method_v1_freezeout_comparison_v2/` and
  `method_v1_freezeout_q0_reference_v3/` provide the ten-point finite-q/q=0
  route comparison and its q=0 reference.
- `fig4_like_freezeout_ratio_dense_20260905_v3/` provides the dense ratio
  comparison, shell tables, audit tables, and figures.
- `freezeout_on_historical_trho_20260905/` provides the historical T-rho
  overlay and its context plots.

The complete main-tree analysis root currently contains 3,582 files and
33,924,685 bytes. Its aggregate inventory hash is recorded in `manifest.json`;
the other retained diagnostic directories are not silently treated as part of
the 32-file light package.

## External worktree decision

`D:/w/jrt-ord` is retained. The four detailed evidence groups match the main
tree byte-for-byte and have the same per-file hashes and aggregate inventory
hash. Its larger root also contains historical/repeated diagnostics that are
not present in main. No worktree or diagnostic file is deleted by this change.

Before deleting `jrt-ord`, create an external archive and independently verify
the full-root inventory in `manifest.json`; the current index alone is not a
replacement archive for the 245+ MB worktree.

## Interpretation boundary

`freezeout_20260907_v2` is an accepted research/smoke-production display
result for the fixed quark-only BQS background. It is not a final-state hadron
yield and not a formal numerical regression baseline. The finite-q/q=0 and
historical-overlay files are review evidence only and must not be used to
override the production route or to claim matched experimental validation.
