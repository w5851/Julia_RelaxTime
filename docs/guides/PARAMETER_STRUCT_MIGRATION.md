# Parameter Struct Migration Guide

> **✅ 迁移状态：已完成（2026-01-25）**
> 参数结构体迁移已全面完成，所有模块支持 `QuarkParams` 和 `ThermoParams` 结构体，同时保持向后兼容。
> 本文档作为参考保留，用于理解迁移过程和双接口模式设计。

## Overview

The PNJL model codebase has been migrated to use `QuarkParams` and `ThermoParams` structs as the standard parameter representation, while maintaining full backward compatibility with existing NamedTuple-based code.

This guide explains:
- The dual interface pattern used throughout the codebase
- When to use structs vs NamedTuples
- How to migrate existing code to use structs
- The internal normalization strategy

## Phase B Status (PNJL workflows/scans)

Phase B extends the migration from RelaxationTime chain to PNJL workflows/scans.

### Completed

1. **Workflow normalization adapter extracted**
     - Added: `src/models/workflows/WorkflowParamAdapters.jl`
     - Shared helpers:
         - `normalize_quark_params`
         - `normalize_thermo_params`
         - `as_legacy_inputs`（已弃用，仅兼容用途）

2. **Workflow modules switched to shared adapter**
     - `src/models/workflows/TransportWorkflow.jl`
     - `src/models/workflows/MesonMassWorkflow.jl`

3. **Structured scan config objects added (non-breaking)**
     - Added: `src/pnjl/scans/ScanConfig.jl`
     - Types:
         - `TmuScanConfig`
         - `TrhoScanConfig`
     - New overloads:
         - `run_tmu_scan(config::TmuScanConfig; kwargs...)`
         - `run_trho_scan(config::TrhoScanConfig; kwargs...)`

### Compatibility guarantees

- Existing kwargs-based scan APIs remain unchanged.
- Existing NamedTuple inputs remain fully supported.
- Config-object style is additive; kwargs still have highest precedence.

## Phase C Status (contract tightening & soft deprecation)

Phase C focuses on making input contracts explicit and testable, while preserving runtime compatibility.

### Completed in current batch

1. **Workflow contract validation hardened**
    - `WorkflowParamAdapters.normalize_quark_params` validates `m/μ` flavor triplets (`u/d/s`) are finite reals.
    - `WorkflowParamAdapters.normalize_thermo_params` validates `T/Φ/Φbar/ξ` are finite reals.
    - Validation errors are unified as `ArgumentError` with field-level messages.

2. **Scan entry input validation added**
    - `run_tmu_scan` / `run_trho_scan` now validate key vector inputs and option symbols at entry.
    - Early contract failures now fail fast with readable `ArgumentError`.

3. **Soft deprecation enabled for workflow NamedTuple path**
    - NamedTuple input for workflow adapter normalization still works.
    - Compatibility path now emits `Base.depwarn` and points to `QuarkParams`/`ThermoParams`.

### Recommended usage in Phase C

- Workflows: prefer `QuarkParams` + `ThermoParams`.
- Scans: prefer config-object entry (`TmuScanConfig` / `TrhoScanConfig`) and use kwargs only for local override.

### Minimal replacement examples

Before:

```julia
q = (m=(u=1.5,d=1.5,s=3.0), μ=(u=0.1,d=0.1,s=0.1))
t = (T=0.15, Φ=0.2, Φbar=0.2, ξ=0.0)
```

After:

```julia
using Main.ParameterTypes: QuarkParams, ThermoParams

q = QuarkParams((m=(u=1.5,d=1.5,s=3.0), μ=(u=0.1,d=0.1,s=0.1)))
t = ThermoParams((T=0.15, Φ=0.2, Φbar=0.2, ξ=0.0))
```

## Phase D Status (hotspot-driven optimization)

Phase D shifts from broad migration to targeted optimization.

### Completed in current batch

1. **Profiling baseline script added**
    - `scripts/dev/profile_paramtypes_hotspots.jl`
    - Covers normalization overhead, representative workflow cost, and request-path transport cost.

2. **First optimization target landed**
    - Module: `src/relaxtime/TransportCoefficients.jl`
    - Path: `TransportRequest` -> transport coefficient entry methods
    - Change: removed `as_namedtuple`-based adapter path in request entry, replaced with direct field-view fast path.

3. **Validation and regression completed**
    - `tests/unit/relaxtime/test_transport_coefficients.jl` fully passed after optimization.
    - Related workflow/scan/contract smoke tests passed in combined regression run.

### Quantitative result (latest local run)

- `transport_coeff(req struct)` / `transport_coeff(explicit nt)` ≈ `0.9178`
- Interpreted as request path being about 8% faster in the measured local baseline.

### Migration recommendation after Phase D

- Continue using **struct-first at public entry**, with compatibility NamedTuple only at explicit boundaries.
- Prefer **hotspot-specific optimization** over full-chain rewrite.

## Phase E Status (RelaxationTime normalization coalescing)

Phase E targets duplicate normalization in the RelaxationTime main path.

### Completed in current batch

1. **Single-boundary normalization in RelaxationTime**
    - Added internal kernel `_compute_average_rates_nt(...)`.
    - `compute_average_rates(...)` and `relaxation_times(...)` now normalize at boundary once, then reuse NamedTuple internal path.
    - Removed duplicate normalization statements in `relaxation_times(...)`.

2. **Compatibility and equivalence checks**
    - Added struct-vs-NamedTuple equivalence test for existing-rates path in `tests/unit/relaxtime/test_relaxation_time.jl`.
    - Existing transport and relaxation test suites remain green.

3. **Profiling extension**
    - Updated `scripts/dev/profile_paramtypes_hotspots.jl` to include lightweight RelaxationTime path comparison.
    - Latest local ratio: `relaxation(struct)/relaxation(nt) ≈ 1.0642`.

### Recommendation after Phase E

- Keep current strategy: one normalization at module boundary + NamedTuple internal kernels.
- Next optimization should focus on nested-chain repeated normalization in `AverageScatteringRate`/`TotalCrossSection` internals, not API-level rewrites.

## Phase F Status (ASR/TCS normalization coalescing)

Phase F lands the next hotspot-driven step in `AverageScatteringRate` (ASR) and `TotalCrossSection` (TCS).

### Completed in current batch

1. **ASR internal kernel split and deep-path cleanup**
    - Added `_precompute_cross_section_nt!`, `_design_w0cdf_s_grid_nt`, `_get_sigma_nt`.
    - Public methods keep dual-interface compatibility, but normalize once at boundary.
    - `_omega_integral_5d` now calls `_get_sigma_nt` directly, reducing repeated boundary normalization in hot loops.

2. **TCS single-boundary normalization path**
    - Added `_total_cross_section_nt` as internal kernel.
    - Added explicit NamedTuple overload for `total_cross_section`.
    - Batch/scan helpers now normalize once and reuse normalized parameters.

3. **Validation**
    - `test_average_scattering_rate.jl` passed with new struct-vs-NamedTuple cache-path equivalence check.
    - `test_relaxation_time.jl` and `test_transport_coefficients.jl` remained green.
    - Profiling script remains runnable after Phase F updates.

### Recommendation after Phase F

- Continue using **single-boundary normalization + internal NamedTuple kernels**.
- For TotalCrossSection profiling, use model-ready fixtures (`A` fields + full K coefficient sets) in dedicated benchmark scripts before judging struct/NamedTuple runtime ratio.

## Phase G Status (TCS model-ready benchmark lane)

Phase G formalizes a dedicated, model-ready profiling path for `TotalCrossSection` (TCS).

### Completed in current batch

1. **Dedicated model-ready profiling script added**
    - Added: `scripts/dev/profile_total_cross_section_model_ready.jl`
    - Script builds fixture with:
        - `A` ensured by `ensure_quark_params_has_A`
        - full effective couplings from `calculate_G_from_A` + `calculate_effective_couplings`

2. **TCS benchmark scope stabilized**
    - Baselines now cover:
        - `total_cross_section`
        - `scan_s_dependence`
        - `calculate_all_total_cross_sections`
    - This avoids false failures caused by incomplete generic fixtures.

3. **Smoke-test guardrail added**
    - Added: `tests/unit/relaxtime/test_total_cross_section_model_ready_fixture_smoke.jl`
    - Verifies single-point, scan, and batch-all TCS APIs on model-ready inputs.

### Quantitative result (latest local run)

- `total_cross_section(:uu_to_uu, s0)` ≈ `0.0168 ms/call`
- `scan_s_dependence(4 points)` ≈ `0.0529 ms/call`
- `calculate_all_total_cross_sections` ≈ `1.0785 ms/call`

### Recommendation after Phase G

- Keep TCS performance tracking in dedicated model-ready scripts.
- Keep generic hotspot script focused on stable, low-dependency comparisons.

## The Dual Interface Pattern

All public functions in the RelaxationTime module chain now accept **both** struct and NamedTuple parameters through Julia's Union types:

```julia
function relaxation_times(
    quark_params::Union{NamedTuple, QuarkParams},
    thermo_params::Union{NamedTuple, ThermoParams},
    K_coeffs::NamedTuple;
    kwargs...
)
```

This means you can call functions with either parameter format, and they will produce identical results.

### How It Works

The dual interface uses an **internal normalization strategy**:

1. **Public API Layer**: Functions accept `Union{NamedTuple, QuarkParams}` and `Union{NamedTuple, ThermoParams}`
2. **Normalization Layer**: Helper functions `_nt_quark(q)` and `_nt_thermo(t)` convert structs to NamedTuples at function entry
3. **Internal Implementation**: All internal logic uses consistent NamedTuple representation

```
┌─────────────────────────────────────────────────────────────┐
│                    Public API Layer                          │
│  Functions accept Union{NamedTuple, QuarkParams/ThermoParams}│
└─────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────┐
│              Normalization Layer                             │
│  _nt_quark(q) / _nt_thermo(t) convert to NamedTuple         │
└─────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────┐
│              Internal Implementation                         │
│  All internal logic uses consistent NamedTuple representation│
└─────────────────────────────────────────────────────────────┘
```

This design ensures:
- **Zero breaking changes**: Existing NamedTuple code continues to work
- **Type stability**: Internal code sees consistent types
- **Zero overhead**: Normalization helpers are inlined by the compiler

## When to Use Structs vs NamedTuples

### Use Structs (Recommended)

**Structs are the recommended pattern for new code** because they provide:

1. **Type safety**: Catch missing fields at construction time
2. **Better IDE support**: Autocomplete and type hints
3. **Clearer intent**: Explicit parameter types in function signatures
4. **Validation**: Constructors can validate field values
5. **Documentation**: Self-documenting parameter structure

**Example:**
```julia
using Main.ParameterTypes: QuarkParams, ThermoParams

# Create struct parameters
quark_params = QuarkParams(
    m=(u=1.52, d=1.52, s=3.04),
    μ=(u=0.3, d=0.3, s=0.3)
)

thermo_params = ThermoParams(0.15, 0.5, 0.5, 0.0)

# Use in calculations
result = relaxation_times(quark_params, thermo_params, K_coeffs; densities=densities)
```

### Use NamedTuples (Backward Compatible)

**NamedTuples are still fully supported** for:

1. **Existing code**: No need to modify working scripts
2. **Quick prototyping**: Faster to write inline
3. **Dynamic construction**: When building parameters programmatically
4. **Legacy interfaces**: When interfacing with older code

**Example:**
```julia
# Create NamedTuple parameters (still works!)
quark_params = (
    m=(u=1.52, d=1.52, s=3.04),
    μ=(u=0.3, d=0.3, s=0.3)
)

thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)

# Use in calculations (identical results)
result = relaxation_times(quark_params, thermo_params, K_coeffs; densities=densities)
```

### Mixed Usage

You can **mix structs and NamedTuples** in the same workflow:

```julia
# Struct for quark params, NamedTuple for thermo params
quark_params = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)

result = relaxation_times(quark_params, thermo_params, K_coeffs; densities=densities)
```

The normalization layer handles all conversions transparently.

## Migrating Existing Code

### Step 1: Import Parameter Types

Add imports at the top of your script:

```julia
using Main.ParameterTypes: QuarkParams, ThermoParams, as_namedtuple
```

### Step 2: Convert NamedTuple Construction to Struct Construction

**Before:**
```julia
quark_params = (
    m=(u=1.52, d=1.52, s=3.04),
    μ=(u=0.3, d=0.3, s=0.3)
)

thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
```

**After:**
```julia
quark_params = QuarkParams(
    m=(u=1.52, d=1.52, s=3.04),
    μ=(u=0.3, d=0.3, s=0.3)
)

thermo_params = ThermoParams(0.15, 0.5, 0.5, 0.0)
```

### Step 3: Update Function Calls (Optional)

Function calls remain **exactly the same**:

```julia
# This works with both structs and NamedTuples
result = relaxation_times(quark_params, thermo_params, K_coeffs; densities=densities)
```

### Step 4: Convert Back to NamedTuple If Needed

If you need to pass parameters to legacy code that only accepts NamedTuples:

```julia
# Convert struct to NamedTuple
quark_nt = as_namedtuple(quark_params)
thermo_nt = as_namedtuple(thermo_params)

# Use with legacy code
legacy_function(quark_nt, thermo_nt)
```

## Complete Migration Example

Here's a complete example showing migration of a typical workflow:

### Before (NamedTuple-based)

```julia
# Old code using NamedTuples
function calculate_transport_properties(T, μ_B)
    # Construct parameters
    quark_params = (
        m=(u=1.52, d=1.52, s=3.04),
        μ=(u=μ_B/3, d=μ_B/3, s=μ_B/3)
    )
    
    thermo_params = (T=T, Φ=0.5, Φbar=0.5, ξ=0.0)
    
    K_coeffs = (
        K_pi=3.67, K_K=5.50, K_eta=0.0,
        K_sigma=3.67, K_zeta=0.0, K_delta=0.0, K_kappa=5.50
    )
    
    # Calculate densities (simplified)
    densities = (
        u=0.1, d=0.1, s=0.05,
        ubar=0.08, dbar=0.08, sbar=0.04
    )
    
    # Compute relaxation times
    result = relaxation_times(quark_params, thermo_params, K_coeffs; densities=densities)
    
    return result
end
```

### After (Struct-based)

```julia
using Main.ParameterTypes: QuarkParams, ThermoParams

# New code using structs
function calculate_transport_properties(T, μ_B)
    # Construct parameters with type safety
    quark_params = QuarkParams(
        m=(u=1.52, d=1.52, s=3.04),
        μ=(u=μ_B/3, d=μ_B/3, s=μ_B/3)
    )
    
    thermo_params = ThermoParams(T, 0.5, 0.5, 0.0)
    
    K_coeffs = (
        K_pi=3.67, K_K=5.50, K_eta=0.0,
        K_sigma=3.67, K_zeta=0.0, K_delta=0.0, K_kappa=5.50
    )
    
    # Calculate densities (simplified)
    densities = (
        u=0.1, d=0.1, s=0.05,
        ubar=0.08, dbar=0.08, sbar=0.04
    )
    
    # Compute relaxation times (same call!)
    result = relaxation_times(quark_params, thermo_params, K_coeffs; densities=densities)
    
    return result
end
```

**Key changes:**
1. Added `using Main.ParameterTypes: QuarkParams, ThermoParams`
2. Changed `quark_params = (...)` to `quark_params = QuarkParams(...)`
3. Changed `thermo_params = (...)` to `thermo_params = ThermoParams(...)`
4. Function call remains identical

## Affected Modules

The following modules support the dual interface:

### Core Modules
- **RelaxationTime.jl**: Main entry point for relaxation time calculations
- **AverageScatteringRate.jl**: Average scattering rate computations
- **TotalCrossSection.jl**: Total cross-section integration
- **ScatteringAmplitude.jl**: Scattering amplitude calculations
- **DifferentialCrossSection.jl**: Differential cross-section calculations
- **TotalPropagator.jl**: Meson propagator calculations

### Utility Modules
- **ParticleSymbols.jl**: Particle symbol parsing and parameter lookup

### Already Supported
- **TransportCoefficients.jl**: Already supports structs via `TransportRequest`

## Internal Normalization Strategy

For developers working on the codebase, here's how the normalization strategy works:

### Normalization Helpers

Each module defines inline normalization helpers:

```julia
@inline _nt_quark(q) = q isa QuarkParams ? as_namedtuple(q) : q
@inline _nt_thermo(t) = t isa ThermoParams ? as_namedtuple(t) : t
```

These helpers:
- Check if input is a struct using `isa`
- Convert to NamedTuple using `as_namedtuple` if needed
- Pass through NamedTuples unchanged
- Are inlined by the compiler for zero overhead

### Function Entry Pattern

All public functions follow this pattern:

```julia
function some_function(
    quark_params::Union{NamedTuple, QuarkParams},
    thermo_params::Union{NamedTuple, ThermoParams},
    ...
)
    # Normalize at entry
    quark_params = _nt_quark(quark_params)
    thermo_params = _nt_thermo(thermo_params)
    
    # Rest of implementation uses NamedTuples
    ...
end
```

This ensures:
- **Type stability**: All downstream code sees consistent NamedTuple types
- **No breaking changes**: Internal implementation unchanged
- **Zero overhead**: Inlining eliminates function call overhead

### Why NamedTuples Internally?

The internal implementation uses NamedTuples (not structs) because:

1. **Backward compatibility**: Existing internal code already uses NamedTuples
2. **Flexibility**: Easy to add/remove fields dynamically (e.g., adding `A` field)
3. **Minimal changes**: No need to rewrite internal logic
4. **Type stability**: Julia's compiler optimizes NamedTuple field access well

The struct interface is purely for the **public API** - it provides type safety and better ergonomics for users while keeping internal implementation simple.

## Performance Considerations

### Zero Overhead Abstraction

The dual interface has **zero runtime overhead** because:

1. **Inlined normalization**: `@inline` ensures helpers are compiled away
2. **Type specialization**: Julia compiles separate versions for struct and NamedTuple inputs
3. **Constant propagation**: Type checks are resolved at compile time

### Benchmarks

Performance tests show struct and NamedTuple versions perform identically:

```julia
using BenchmarkTools

q_struct = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
q_nt = as_namedtuple(q_struct)

# Both versions have identical performance
@btime relaxation_times($q_struct, $t_struct, $K; densities=$densities)
@btime relaxation_times($q_nt, $t_nt, $K; densities=$densities)
```

## Error Handling

### Struct Construction Errors

Structs provide early error detection:

```julia
# Missing required field
try
    q = QuarkParams((m=(u=1.52, d=1.52, s=3.04)))  # Missing μ
catch e
    # Error: "QuarkParams: input is missing field :μ"
end

# Missing thermodynamic parameter
try
    t = ThermoParams((T=0.15, Φ=0.5))  # Missing Φbar
catch e
    # Error: "ThermoParams: input is missing field :Φbar"
end
```

### NamedTuple Field Access Errors

NamedTuples detect errors at field access time:

```julia
# Missing field in NamedTuple
quark_params = (m=(u=1.52, d=1.52, s=3.04),)  # Missing μ

try
    result = relaxation_times(quark_params, thermo_params, K_coeffs; densities=densities)
catch e
    # Error: "type NamedTuple has no field μ"
end
```

**Recommendation**: Use structs to catch errors earlier (at construction time rather than usage time).

## FAQ

### Q: Do I need to update my existing scripts?

**A:** No! All existing NamedTuple-based code continues to work without modification. The migration is fully backward compatible.

### Q: What's the benefit of using structs?

**A:** Structs provide type safety, better IDE support, earlier error detection, and clearer documentation. They're the recommended pattern for new code.

### Q: Can I mix structs and NamedTuples?

**A:** Yes! You can use structs for some parameters and NamedTuples for others in the same function call. The normalization layer handles all conversions.

### Q: Is there a performance difference?

**A:** No. The dual interface has zero runtime overhead due to inlining and type specialization. Benchmarks show identical performance.

### Q: How do I convert between structs and NamedTuples?

**A:** Use `as_namedtuple(struct)` to convert structs to NamedTuples. To convert back, use the struct constructor: `QuarkParams(namedtuple)`.

### Q: What if I need to add the `A` field?

**A:** Use `ensure_quark_params_has_A(quark_params, thermo_params)`. This function works with both structs and NamedTuples and returns a NamedTuple with the `A` field added.

### Q: Which modules support the dual interface?

**A:** All modules in the RelaxationTime chain: RelaxationTime, AverageScatteringRate, TotalCrossSection, ScatteringAmplitude, DifferentialCrossSection, TotalPropagator, and ParticleSymbols.

## Summary

The parameter struct migration provides:

✅ **Type safety** through struct constructors  
✅ **Backward compatibility** with existing NamedTuple code  
✅ **Zero overhead** through inlined normalization  
✅ **Flexible usage** with mixed struct/NamedTuple support  
✅ **Better documentation** through explicit types  

**Recommendation**: Use structs for new code, but don't worry about updating existing code unless you want the benefits of type safety.
