# Other-Mott local slope diagnostics

This case tests local transport slope modulation near other meson Mott thresholds, using the frozen v15 input grid. The planned finite-window diagnostic review is complete as of 2026-10-31. K windows also show slope modulation, including windows without overlapping pi brackets. The narrowest pi slopes and ambiguous mixed-meson branches retain explicit limits. [The Mott analysis](../../../../../../../docs/analysis/relaxtime/threshold_cache_upgrade_20261021_analysis/mott_local_bumps_20261028.md) maintains the interpretation; [scientific_review_v1.json](scientific_review_v1.json) binds this review to the returned evidence. No formal data, publication pointer or PDF is changed by this case.

| Runtime field | Value |
| --- | --- |
| Compute host | `quote-remote` |
| SSH configuration | `D:/Desktop/CodexBountyLab/executions/cross-chain-quote-observer-2026-10-15/tunnel-setup-2026-10-16/ssh_config` |
| Remote case | `D:/CodexCompute/trmott-other-slopes-20261029-v1` |
| Frozen numerical source | `D:/CodexCompute/trfull-20261026-v1/case/source` |
| Numerical commit | `2e5c411c716d5cbc524c3fc6cfb335fe05b39444` |
| Completed recovery task | `CodexMottOtherSlopes-20261030-recovery-v1` |
| Completed reference audit task | `CodexMottOtherSlopes-20261030-final-acceptance-v1` |
| Queue states | Remote `recovery_v1.status.json` and `final_acceptance_v1.status.json` |
| Frozen queue configuration | [pipeline_v2.config.json](pipeline_setup_v2/pipeline_v2.config.json) |

The [complete reference review](complete_acceptance_snapshot_v1/returned/acceptance_review_complete_v1/summary.json) validates 195 main points, 3315 rates and all 175 planned reference evaluations: 127 original-coordinate controls and 48 new xi points. It includes 819 matched-endpoint slope comparisons and 432 xi-split comparisons without missing intervals. The [733-file return](complete_acceptance_snapshot_v1/local_return_verification.json) and archive hashes were checked again during continuation. The remote recovery and final-audit states remain `ready_for_agent_review`, with no active, pending or failed runs; both finished at 13:34 Beijing time on 2026-10-31.

The [same-width decomposition](center_modulation_acceptance_v1/returned/summary.json) compares actual symmetric triplets at identical half-width h, preserving endpoints when reference midpoints are present. It verifies 1035 centered modulations against the frozen review and separates process-rate/density and tau/final-coefficient contributions. The maximum decomposition closure error is 5.95790e-11. Its [nine-file return](center_modulation_acceptance_v1/local_return_verification.json), source snapshots and 33 input hash bindings are verified; [interpretation_views.json](center_modulation_acceptance_v1/returned/interpretation_views.json) contains the remotely computed window and species tables. No interpolation or new physics solve was introduced for this decomposition.

The [complete PNG index](../../phase_guided/mott_other_slopes_20261031_v1_complete_png_review/README.md) contains 65 color and 65 grayscale plots for tau_ubar, tau_sbar, eta/s, zeta/s and sigma/T. The manifest contract and all output hashes pass. [The visual review](../../phase_guided/mott_other_slopes_20261031_v1_complete_png_review/agent_review.json) distinguishes auxiliary contact-sheet screening, primary-agent detail checks and machine-only coverage. The auxiliary model supplies layout candidates only; scientific conclusions use remote tables and source. All plots retain the original values and secants, require 6.75-inch width, and remain diagnostic PNGs awaiting author review.

Acceptance is limited to the planned finite xi windows under the inherited integration settings. Tighter s tables do not certify all integration layers or a continuous derivative limit. The pi controls still change substantially between h=0.0002 and h=0.0001. The K process groups retain other exchanges and interference; overlapping pi brackets in other02/04/13 preclude exclusive attribution. The sigma window contains background trends and weaker central details, not evidence of exactly zero response. Existing eta/eta-prime root ambiguities remain unresolved. The frozen automatic `whole_task_complete=false` fields are preserved; the separate scientific review records the completed scope without promoting data or figures. v15 PDF remains deferred by the author.

## Retained earlier stages and recovery

The [618-file main-grid snapshot](main_grid_snapshot_v1/local_return_verification.json) retains its historical 43/175-reference scope, 198 matched-endpoint comparisons and 108 xi splits. Its [39-pair PNG manifest](../../../../../figures/relaxtime/transport/phase_guided/mott_other_slopes_20261030_v1_main_grid_png_review/plot_manifest.json) and [visual review](../../phase_guided/mott_other_slopes_20261030_v1_main_grid_png_review/agent_review.json) remain unchanged. They do not describe the final reference coverage.

The original v2 queue stopped dispatching after the first other05 tighter point exceeded its 80000-evaluation cap; already running workers finished naturally. Recovery reused completed runs and raised only max_evaluations to 320000 and max_nodes to 96000, retaining the 1e-6 tolerance and physical quadratures. The failed point passed, with its two expensive channels requiring 124274 and 108096 evaluations. The [configuration and failure evidence](main_grid_snapshot_v1/returned/acceptance_evidence_v1/) remain available. Five recovery and three acceptance unit tests passed on quote-remote. The earlier v1 controller and its [supersession receipt](pipeline_setup_v2/pipeline_v1_superseded.json) also remain available.

The separate final task waited for an idle successful recovery state before running the frozen acceptance script with `--require-complete`; it did not start numerical workers. The [frozen setup](final_acceptance_setup_v1/config.json), [source/hash preflight](final_acceptance_setup_v1/preflight.json) and [five guard tests](final_acceptance_setup_v1/finalization_unit_v1.stderr.log) retain its execution contract. Do not start another queue copy when resuming. A read-only status check is:

```powershell
ssh -F 'D:/Desktop/CodexBountyLab/executions/cross-chain-quote-observer-2026-10-15/tunnel-setup-2026-10-16/ssh_config' quote-remote powershell -NoProfile -NonInteractive -Command "Get-Content -LiteralPath 'D:/CodexCompute/trmott-other-slopes-20261029-v1/recovery_v1.status.json' -Raw"
```

The complete package is already downloaded and verified; do not extract over its existing `returned/` directory. If a fresh copy is required, use a new destination and `scripts/analysis/relaxtime/verify_transport_other_mott_return.py`, which verifies/extracts bytes and metadata only. All physics calculations, derived numerical analyses and numerical tests run on quote-remote. Any future pi-center convergence work must use a new case and retain the present results.

The [449-file pilot precision return](pilot_precision_snapshot_v1/local_return_verification.json) retains its four-point tighter pilot, 27 matched-endpoint slope comparisons and selection plan. Its [five-pair PNG manifest](../../../../../figures/relaxtime/transport/phase_guided/mott_other_slopes_20261029_v1_pilot_precision_png_review/plot_manifest.json), the earlier [410-file pilot return](pilot_snapshot_v1/local_return_verification.json) and [six-point pilot PNG manifest](../../../../../figures/relaxtime/transport/phase_guided/mott_other_slopes_20261029_v1_pilot_png_review/plot_manifest.json) retain their historical scope and bytes.
