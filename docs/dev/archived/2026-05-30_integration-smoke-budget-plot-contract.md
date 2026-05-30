---
title: Integration Smoke Budget: Full-Only Demotions
archived: true
original: docs/dev/active/2026-05-30_integration-smoke-budget-plot-contract.md
archived_date: 2026-05-30
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Integration Smoke Budget: Full-Only Demotions

## Context

`tests/integration/runtests.jl` defines the integration smoke profile as an ultra-fast local edit-run subset with a warm-cache target near one minute. Several existing smoke files launch Julia subprocesses, run orchestrator paths, or perform heavier scan setup, so they no longer fit the smoke budget.

## Scope

- Demote local-overbudget files out of `INTEGRATION_SMOKE_FILES`.
- Keep the test covered by `INTEGRATION_PROFILE=full`, where all integration `test_*.jl` files are included.
- Preserve the remaining relaxtime smoke coverage in the smoke list.

## Acceptance

- [x] Smoke list no longer includes `test_plot_contract_smoke.jl`.
- [x] Smoke list no longer includes other measured local-overbudget subprocess/orchestrator files.
- [x] `INTEGRATION_PROFILE=smoke` completes within the local smoke budget.
- [x] Demoted files remain covered structurally by `INTEGRATION_PROFILE=full` directory scanning.

## Timing Evidence

- `config/test_config_profile_smoke.jl`: 34.6s.
- `relaxtime/test_plot_contract_smoke.jl`: did not finish within 120s.
- `relaxtime/test_cross_section_orchestrated_smoke.jl`: 156.4s.
- `relaxtime/test_xi_smoothness_batch_runner_smoke.jl`: 33.1s.
- `relaxtime/test_phase_guided_transport_scan_smoke.jl`: 66.1s.
- After demotion, `INTEGRATION_PROFILE=smoke`: 18 pass, 16.6s.
