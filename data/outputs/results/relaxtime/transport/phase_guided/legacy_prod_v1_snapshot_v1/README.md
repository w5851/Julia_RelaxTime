# RS `prod_v1` legacy result snapshot

This versioned directory is the byte-preserving snapshot of the superseded
Issue #130 RS transport result trees.  The old trees no longer occupy the
mode-specific canonical roots; no numerical file was recomputed or deleted.

- canonical runtime result: `../../<mode>/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2/`
- snapshot case: `<mode>/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1/`
- use: explicit legacy fallback, rollback, and historical reproduction only
- calculation SHA: `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- source commit: `05be2c05186f8e12baf3097b68f8619e53d19711`

`RETIREMENT_MANIFEST.json` records the original canonical path, snapshot path,
per-tree byte count, inner manifest hash, and the path-independent tree hash.
The inner case manifests and convergence sidecars are intentionally unchanged;
their historical path strings are provenance, not active runtime locators.

Physical deletion is not authorized by this snapshot migration and requires a
separate review and PR.
