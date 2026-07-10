---
name: api-doc-authoring
description: "Author or restructure Julia_RelaxTime public API documentation under docs/api using three views: user entrypoints, responsibility core, and an automatically generated complete export index. Use for stable public APIs and export coverage; do not use for general design documents, paper prose, or implementation-task execution."
---

# API Documentation Authoring

## Core contract

- Organize documentation by user need, not source-file order.
- Put stable façades and complete workflows before internal numerical helpers.
- Maintain user entrypoints and responsibility-core explanations manually.
- Generate the complete export index with `scripts/dev/generate_api_export_index.jl`; never maintain it by hand.
- Update `docs/api/` whenever a stable public entrypoint changes.
- Preserve units, stability labels, compatibility notes, and repository terminology.

## Evidence sources

Inspect the target export surface and its authoritative callers before drafting:

- `src/models/Models.jl` and `src/models/entrypoints.jl`;
- target module `export` lists and aggregation modules;
- `README.md`, existing `docs/api/`, guides, scripts, HTTP routes, and tests.

Treat an exported symbol as a public candidate, not automatically as a recommended user entrypoint.

## Workflow

1. Build a deduplicated export inventory from the target module and its aggregators.
2. Identify stable user entrypoints from the unified façade, README, scripts, server routes, and existing API docs.
3. Identify responsibility-core material: algorithms, criteria, data flow, and module boundaries needed by maintainers.
4. Choose a single-page layout for small topics or a topic directory for larger surfaces.
5. Draft the user view with positioning, units, input contract, stability, and the shortest runnable example.
6. Draft the responsibility-core view without copying implementation files in source order.
7. Run the export-index generator and treat its output as the completeness baseline.
8. Mark exported-but-discouraged or experimental symbols explicitly instead of omitting them.
9. Verify links, symbol coverage, terminology, units, and migration notes.

## Classification heuristics

A symbol is likely a user entrypoint when it is re-exported by a top-level aggregator, documented in quickstarts, called by stable scripts or HTTP routes, and completes a workflow rather than a local calculation.

A topic belongs in the responsibility core when it implements a key physical/numerical criterion, controls branching or data flow, or is reused by multiple public workflows.

When confidence is low, use the conservative default: promote only entrypoints explicitly recommended by the unified façade or README; place the rest in the generated export index until evidence supports a stronger classification.

## Output structure

Include:

- module purpose and stability;
- units, inputs, outputs, and dependencies;
- recommended user entrypoints and examples;
- responsibility-core concepts and boundaries;
- generated complete export index;
- compatibility and migration notes.

For `Models` topic reorganization, old PNJL/relaxtime page migration, magnetic variants, derived quantities, or transport placement, read [references/models-taxonomy.md](references/models-taxonomy.md).

## Validation

- Run `scripts/dev/generate_api_export_index.jl` for the target surface.
- Run `scripts/dev/check_docs_consistency.jl` and relevant entry-contract checks.
- Confirm the generated index contains every export and records whether human-authored docs mention it.
- Confirm examples use current `Models` entrypoints rather than removed implementation paths.
