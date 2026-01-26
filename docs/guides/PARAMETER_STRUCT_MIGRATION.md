# Parameter Struct Migration Guide

## Overview

The PNJL model codebase has been migrated to use `QuarkParams` and `ThermoParams` structs as the standard parameter representation, while maintaining full backward compatibility with existing NamedTuple-based code.

This guide explains:
- The dual interface pattern used throughout the codebase
- When to use structs vs NamedTuples
- How to migrate existing code to use structs
- The internal normalization strategy

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
