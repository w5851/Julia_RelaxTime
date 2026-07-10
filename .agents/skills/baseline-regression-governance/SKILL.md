---
name: baseline-regression-governance
description: "Manage repository-wide numerical baseline governance: baseline generation, storage, regression assertions, CI tiers, tolerance evidence, and baseline-change admission. Use for cross-domain policy or new baseline systems; use transport-regression-keeper instead for day-to-day transport/relaxtime drift diagnosis."
---

# Baseline Regression Governance

## Invariants

- Treat a baseline as a versioned numerical contract, not a convenient output snapshot.
- Never refresh a baseline merely to make a failing test pass.
- Explain every tolerance through numerical error, physical meaning, or validation evidence.
- Keep baseline data deterministic and free of timestamps, machine paths, and random ordering.
- Preserve the repository's `Models` entrypoint and solver semantics unless the task explicitly changes them.

## Test-layer decision

- Put fixed-point or sampled numerical promises in `tests/regression/<domain>/`.
- Use `tests/unit/<domain>/` only for baseline readers, serializers, validators, or local numerical helpers.
- Use `tests/integration/<domain>/` for cross-module workflow assembly without a long-lived numerical target.
- Use `tests/validation/<domain>/` for external references, literature mappings, or acceptance criteria.
- Keep a small deterministic regression set in smoke; place broad grids or expensive matrices in core/full/nightly coverage.

## Repository layout

- Baselines: `tests/baselines/<domain>/baseline_<feature>_<scope>_vN.csv`
- Regression assertions: `tests/regression/<domain>/test_<feature>_regression.jl`
- Baseline export or comparison scripts: `scripts/dev/export_<feature>_baseline.jl`

Keep point identifiers, compared observables, units, and column order explicit. Use stable numeric formatting and avoid volatile metadata in CSV files.

## Workflow

1. Define the protected behavior, fixed points, observables, units, and expected precision.
2. Select the test layer with the decision rules above.
3. Generate candidate values through a reproducible script and record the source commit/configuration.
4. Compare the candidate with the current baseline point by point.
5. Diagnose drift before deciding whether code, expectations, or the baseline is wrong.
6. Add or update the narrowest sufficient regression coverage.
7. Run the focused selector first, then expand to the appropriate profile and governance checks.
8. Record the reason, impact, numerical delta, commands, and residual risk for every accepted baseline change.

## Change admission

An accepted baseline change must include:

- the semantic or numerical reason;
- affected modules, points, and observables;
- an old/new difference summary;
- unchanged constraints that were rechecked;
- validation commands and results;
- relevant documentation or active-task updates.

Reject unexplained drift, tolerance-only fixes, and baseline-only commits without evidence.

## Current transport example

- Baseline: `tests/baselines/relaxtime/baseline_transport_fixedpoints_v1.csv`
- Export script: `scripts/dev/export_transport_fixedpoint_baseline.jl`
- Regression test: `tests/regression/relaxtime/test_transport_fixedpoint_regression.jl`
- Focused run:

```powershell
julia --project=. -e 'ENV["REGRESSION_FILES"]="relaxtime/test_transport_fixedpoint_regression.jl"; include("tests/regression/runtests.jl")'
```

Use `transport-regression-keeper` to interpret transport-specific drift and choose any additional integration or validation coverage.
