---
name: julia-pro
description: Julia 1.10+ engineering guidance for scientific/numerical code, performance optimization, multiple dispatch design, and package-quality implementation. Use when tasks involve Julia code structure, type stability, profiling, testing, or refactoring in computational backends.
license: MIT
metadata:
  author: local-adapted
  version: "1.0.0"
---

# Julia Pro

## When to apply
- Implementing or refactoring Julia scientific backend code.
- Diagnosing performance issues (allocations, type instability, slow hot paths).
- Designing extensible APIs with multiple dispatch.
- Improving package quality (tests, docs, reproducibility).

## Core workflow
1. Clarify numerical goal, invariants, and expected precision/tolerance.
2. Make data contracts explicit at boundaries (input/output types and units).
3. Keep kernels type-stable and allocation-light.
4. Validate correctness before optimization.
5. Profile hot paths, then optimize only measured bottlenecks.

## Implementation rules
- Prefer concrete container element types in hot loops.
- Avoid global mutable state in compute kernels.
- Isolate IO/config parsing from numerical kernels.
- Use function barriers when dealing with abstract inputs.
- Avoid type piracy; extend methods in your own namespace.
- Keep APIs composable and predictable.

## Performance checklist
- Run allocation checks for representative workloads.
- Confirm return type stability for key kernels.
- Benchmark before/after each optimization change.
- Guard against regressions with deterministic smoke cases.

## Testing expectations
- Add/maintain unit tests for numerical correctness.
- Include edge-case tests for domain boundaries.
- Use tolerance-aware assertions for floating-point outputs.
- Keep fixtures minimal and reproducible.

## Response format
- Summary of approach
- Proposed/implemented code changes
- Verification method and benchmark notes
- Remaining risks and next actions
