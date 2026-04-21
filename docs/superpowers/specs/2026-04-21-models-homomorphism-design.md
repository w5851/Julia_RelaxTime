# Models Homomorphism Design (#96)

## Context

- #92 is closed after boundary convergence, export layering hardening, and workflow naming unification.
- New target (#96): make physical-model structure in `src/models` homomorphic across all model kinds, and align API contracts in the same round.
- User constraints for #96:
  - Do directory + API together.
  - Full-model one-shot scope (`njl`, `njl2`, `pnjl`, `rpnjl`, `pnjl_magnetic`, `rotation`, `gas_liquid`).
  - Remove old-path compatibility includes in this round.

## Goals

1. Unify per-model directory semantics.
2. Unify minimum model API contract and capability signaling.
3. Remove old include-path compatibility layer in scripts/tests (switch to new canonical paths/API access).
4. Keep `Models` public behavior stable for existing entrypoints.

## Non-Goals

- No physics-formula redefinition in this issue.
- No numeric baseline target changes.
- No broad solver algorithm redesign.

## Canonical Model Skeleton

Each model kind should conform to this semantic skeleton under `src/models/<model_kind>/`:

- `api.jl`: stable model-facing API surface for this model kind.
- `core/`: model core physics/numerical implementation.
- `adapters/`: mapping to unified solver/workflow/transport contracts.
- `workflows/`: model-specific workflows (optional by capability, but directory semantics retained).
- `capabilities.jl`: explicit capability declaration and missing-capability signaling policy.

For model kinds that already have deeper internals (for example `pnjl_physics/core`), migration keeps internals but introduces/normalizes top-level semantic anchors above.

## Unified Contract

Minimum contract symbols to unify by model kind capability:

- `solve_gap(model, T, mu; kwargs...)`
- `model_thermo(model, x_state, mu_vec, T; kwargs...)`
- `number_densities(model, x_state, T, mu_vec; kwargs...)`

Capability policy:

- Introduce/normalize `UnsupportedCapabilityError` usage where capability is absent.
- Missing support must fail through explicit capability checks, not through missing method accidents.

## Execution Strategy (Approved: Plan A)

Although scope is one-shot, execution remains internally staged for control:

### Stage 1: Structure Homomorphism

- Build canonical per-model anchors (`api/core/adapters/workflows/capabilities`) for all 7 model kinds.
- Move/repoint files to canonical locations.
- Update `Models.jl` includes to canonical paths only.

### Stage 2: API Homomorphism

- Normalize per-model API entry wrappers in `api.jl`.
- Ensure minimum contract is present (or explicitly unsupported by capability rule).
- Centralize capability declaration and contract checks.

### Stage 3: Repo-Wide Cutover

- Replace scripts/tests old include paths and direct legacy model-file path assumptions.
- Remove old compatibility include layer in this round (per user requirement).
- Re-run governance and layered tests.

## Verification Plan

### Governance (must pass)

- `julia --project=. scripts/dev/check_models_entry_contract.jl`
- `julia --project=. scripts/dev/check_solver_contract_leakage.jl`
- `julia --project=. scripts/dev/check_pnjl_migration_guard.jl`

### Test gates (minimum)

- Targeted unit/integration around model entrypoints and workflows.
- Regression/validation smoke gates covering renamed/moved model paths.

### Acceptance Evidence

- Before/after model structure matrix in #96 comments.
- Capability matrix in #96 comments.
- CI evidence links for final one-shot batch.

## Risks and Mitigations

1. **Include path breakage in scripts/tests**
   - Mitigation: exhaustive path grep + targeted smoke runs before full CI.
2. **Implicit capability assumptions**
   - Mitigation: explicit capability declarations and deterministic `UnsupportedCapabilityError` path.
3. **One-shot blast radius**
   - Mitigation: keep internal staged validation checkpoints, but single thematic delivery.

## Rollback Strategy

- Keep changes in coherent commit slices by stage.
- If Stage 3 regressions occur, revert latest slice only, preserve Stage 1/2 where safe.

## Definition of Done for #96

- All 7 model kinds conform to canonical directory semantics.
- Minimum API contract is unified or explicitly unsupported via capability policy.
- Old compatibility include paths removed for model structure layer.
- Governance scripts pass.
- Relevant layered tests pass.
