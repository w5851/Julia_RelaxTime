# Parameter Types API Documentation

This document provides comprehensive API documentation for the `QuarkParams` and `ThermoParams` structs and their associated helper functions.

## Table of Contents

- [QuarkParams](#quarkparams)
- [ThermoParams](#thermoparams)
- [Helper Functions](#helper-functions)
- [Normalization Helpers](#normalization-helpers)
- [Phase B Additions](#phase-b-additions)
- [Phase C Additions](#phase-c-additions)
- [Phase D Additions](#phase-d-additions)
- [Phase E Additions](#phase-e-additions)
- [Phase F Additions](#phase-f-additions)
- [Phase G Additions](#phase-g-additions)
- [Usage Examples](#usage-examples)

---

## QuarkParams

### Description

`QuarkParams` is a struct that encapsulates quark masses and chemical potentials for the three quark flavors (u, d, s) in the PNJL model.

### Definition

```julia
struct QuarkParams
    m::NamedTuple  # (u=..., d=..., s=...)
    μ::NamedTuple  # (u=..., d=..., s=...)
end
```

### Fields

- **`m`**: NamedTuple containing quark masses in fm⁻¹
  - `m.u`: Up quark mass
  - `m.d`: Down quark mass
  - `m.s`: Strange quark mass

- **`μ`**: NamedTuple containing quark chemical potentials in fm⁻¹
  - `μ.u`: Up quark chemical potential
  - `μ.d`: Down quark chemical potential
  - `μ.s`: Strange quark chemical potential

### Construction

#### From NamedTuple

```julia
q = QuarkParams((
    m=(u=1.52, d=1.52, s=3.04),
    μ=(u=0.3, d=0.3, s=0.3)
))
```

#### From Keyword Arguments (if constructor added)

```julia
q = QuarkParams(
    m=(u=1.52, d=1.52, s=3.04),
    μ=(u=0.3, d=0.3, s=0.3)
)
```

### Validation

The constructor validates that both `m` and `μ` fields are present. Missing fields will raise an error:

```julia
# This will fail
try
    q = QuarkParams((m=(u=1.52, d=1.52, s=3.04),))  # Missing μ
catch e
    println(e)  # "QuarkParams: input is missing field :μ"
end
```

### Field Access

```julia
q = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))

# Access masses
m_u = q.m.u  # 1.52
m_d = q.m.d  # 1.52
m_s = q.m.s  # 3.04

# Access chemical potentials
μ_u = q.μ.u  # 0.3
μ_d = q.μ.d  # 0.3
μ_s = q.μ.s  # 0.3
```

---

## ThermoParams

### Description

`ThermoParams` is a struct that encapsulates thermodynamic parameters for the PNJL model, including temperature, Polyakov loops, and anisotropy.

### Definition

```julia
struct ThermoParams
    T::Float64      # Temperature
    Φ::Float64      # Polyakov loop
    Φbar::Float64   # Conjugate Polyakov loop
    ξ::Float64      # Anisotropy parameter
end
```

### Fields

- **`T`**: Temperature in fm⁻¹
- **`Φ`**: Polyakov loop (dimensionless, typically 0-1)
- **`Φbar`**: Conjugate Polyakov loop (dimensionless, typically 0-1)
- **`ξ`**: Anisotropy parameter (dimensionless, 0 for isotropic)

### Construction

#### Direct Construction

```julia
t = ThermoParams(0.15, 0.5, 0.5, 0.0)
```

#### From NamedTuple

```julia
t = ThermoParams((T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0))
```

### Validation

The constructor validates that all four fields are present. Missing fields will raise an error:

```julia
# This will fail
try
    t = ThermoParams((T=0.15, Φ=0.5))  # Missing Φbar and ξ
catch e
    println(e)  # "ThermoParams: input is missing field :Φbar"
end
```

### Field Access

```julia
t = ThermoParams(0.15, 0.5, 0.5, 0.0)

# Access fields
T = t.T      # 0.15
Φ = t.Φ      # 0.5
Φbar = t.Φbar  # 0.5
ξ = t.ξ      # 0.0
```

---

## Helper Functions

### `as_namedtuple`

Convert a struct to its NamedTuple representation.

#### Signature

```julia
as_namedtuple(q::QuarkParams)::NamedTuple
as_namedtuple(t::ThermoParams)::NamedTuple
```

#### Arguments

- `q`: QuarkParams struct to convert
- `t`: ThermoParams struct to convert

#### Returns

A NamedTuple with the same fields as the input struct.

#### Examples

```julia
# Convert QuarkParams to NamedTuple
q = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
q_nt = as_namedtuple(q)
# q_nt = (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))

# Convert ThermoParams to NamedTuple
t = ThermoParams(0.15, 0.5, 0.5, 0.0)
t_nt = as_namedtuple(t)
# t_nt = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
```

#### Use Cases

- Interfacing with legacy code that expects NamedTuples
- Serialization and deserialization
- Dynamic field manipulation (merge, modify, etc.)

---

## Normalization Helpers

These are internal helper functions used by modules to support the dual interface pattern. They are not exported but are documented here for developers working on the codebase.

### `_nt_quark`

Convert QuarkParams to NamedTuple if needed, otherwise pass through.

#### Signature

```julia
@inline _nt_quark(q) = q isa QuarkParams ? as_namedtuple(q) : q
```

#### Arguments

- `q`: Either a QuarkParams struct or a NamedTuple

#### Returns

A NamedTuple representation of the quark parameters.

#### Behavior

- If `q` is a `QuarkParams` struct, converts it to NamedTuple using `as_namedtuple`
- If `q` is already a NamedTuple, returns it unchanged
- Inlined for zero runtime overhead

#### Usage

```julia
function some_function(quark_params::Union{NamedTuple, QuarkParams})
    # Normalize at function entry
    quark_params = _nt_quark(quark_params)
    
    # Rest of implementation uses NamedTuple
    m_u = quark_params.m.u
    ...
end
```

### `_nt_thermo`

Convert ThermoParams to NamedTuple if needed, otherwise pass through.

#### Signature

```julia
@inline _nt_thermo(t) = t isa ThermoParams ? as_namedtuple(t) : t
```

#### Arguments

- `t`: Either a ThermoParams struct or a NamedTuple

#### Returns

A NamedTuple representation of the thermodynamic parameters.

#### Behavior

- If `t` is a `ThermoParams` struct, converts it to NamedTuple using `as_namedtuple`
- If `t` is already a NamedTuple, returns it unchanged
- Inlined for zero runtime overhead

#### Usage

```julia
function some_function(thermo_params::Union{NamedTuple, ThermoParams})
    # Normalize at function entry
    thermo_params = _nt_thermo(thermo_params)
    
    # Rest of implementation uses NamedTuple
    T = thermo_params.T
    ...
end
```

---

## Phase C Additions

Phase C starts tightening input contracts for PNJL workflows/scans while keeping backward compatibility.

### Recommended entry signatures

- Scans (recommended):
    - `run_tmu_scan(config::TmuScanConfig; kwargs...)`
    - `run_trho_scan(config::TrhoScanConfig; kwargs...)`
- Workflows:
    - `solve_gap_and_transport(T_fm::Real, mu_fm::Real; ...)`
    - `solve_transport_from_equilibrium(base, T_fm::Real, mu_fm::Real; ...)`
    - `solve_gap_and_meson_point(T_fm::Real, mu_fm::Real; ...)`

### Contract validation behavior

- Workflow adapter layer (`WorkflowParamAdapters`) now validates:
    - `quark_params.m.(u,d,s)` and `quark_params.μ.(u,d,s)` are finite real values.
    - `thermo_params.(T,Φ,Φbar,ξ)` are finite real values.
- Scan entry layer validates:
    - `T_values/mu_values/rho_values/xi_values` are non-empty `AbstractVector{<:Real}`.
    - Enum-like options (`seed_policy`, `constraint_mode`, `thermo_backend`, `solver_backend`) are within allowed symbols.
- Invalid contract inputs throw `ArgumentError` with field-level messages.

### Deprecation path (soft migration)

- Passing raw `NamedTuple` to workflow adapter normalization remains supported but now emits `Base.depwarn`.
- Preferred replacement:
    - `NamedTuple -> QuarkParams`
    - `NamedTuple -> ThermoParams`

Minimal migration snippet:

```julia
using Main.ParameterTypes: QuarkParams, ThermoParams

qp = QuarkParams((m=(u=1.5, d=1.5, s=3.0), μ=(u=0.1, d=0.1, s=0.1)))
tp = ThermoParams((T=0.15, Φ=0.2, Φbar=0.2, ξ=0.0))
```

---

## Phase D Additions

Phase D introduces hotspot-driven optimization on top of Phase C contracts.

### Implemented optimization target

- Module: `src/relaxtime/TransportCoefficients.jl`
- Target path: `TransportRequest` entry to
    - `shear_viscosity(req)`
    - `electric_conductivity(req)`
    - `bulk_viscosity_isentropic(req)`
    - `transport_coefficients(req)`

### What changed

- Removed request-entry reliance on `as_namedtuple(req.quark/req.thermo)`.
- Added explicit request fast path using direct field views at call boundary.
- Preserved public compatibility:
    - NamedTuple public APIs remain unchanged.
    - Direct `QuarkParams` + `ThermoParams` overloads remain available.

### Validation summary

- Existing `TransportCoefficients` unit tests pass.
- Added/extended struct-vs-legacy equivalence assertions in tests.
- Phase D profiling script reports request-path improvement in latest run:
    - `transport_coeff(req struct) / transport_coeff(explicit nt) ≈ 0.9178`
    - Interpreted as ~8% faster in this local baseline.

---

## Phase E Additions

Phase E continues hotspot-driven optimization in the RelaxationTime chain.

### Implemented target

- Module: `src/relaxtime/RelaxationTime.jl`
- Target path: `relaxation_times(...)` calling `compute_average_rates(...)`

### What changed

- Added internal NamedTuple kernel: `_compute_average_rates_nt(...)`.
- `compute_average_rates(...)` now normalizes once at boundary and delegates to internal kernel.
- `relaxation_times(...)` now normalizes once and directly reuses the internal kernel, removing repeated normalization calls.

### Validation summary

- New equivalence test added in `test_relaxation_time.jl` for struct vs NamedTuple on existing-rates path.
- Existing relaxation and transport coefficient tests pass.
- Profiling script now includes lightweight `relaxation_times(...; existing_rates=...)` comparison:
    - latest local run ratio `relaxation(struct)/relaxation(nt) ≈ 1.0642`
    - indicates functional parity and that remaining overhead is small, guiding next hotspot selection.

---

## Phase F Additions

Phase F coalesces avoidable normalization in `AverageScatteringRate` and `TotalCrossSection`.

### Implemented changes

- `AverageScatteringRate.jl`
    - Added internal NamedTuple kernels:
        - `_precompute_cross_section_nt!`
        - `_design_w0cdf_s_grid_nt`
        - `_get_sigma_nt`
    - Public wrappers still accept both struct and NamedTuple, but now normalize once at boundary.
    - Deep-loop path (`_omega_integral_5d`) now calls `_get_sigma_nt` directly to avoid repeated normalization calls.

- `TotalCrossSection.jl`
    - Added internal kernel `_total_cross_section_nt`.
    - Added explicit NamedTuple overload for `total_cross_section` to bypass wrapper normalization overhead in internal/batch calls.
    - `calculate_all_total_cross_sections` and `scan_s_dependence` now normalize once and reuse normalized parameters.

### Validation summary

- `test_average_scattering_rate.jl` passed, including new struct-vs-NamedTuple equivalence case with cached sigma path.
- `test_relaxation_time.jl` and `test_transport_coefficients.jl` remained green after refactor.
- Hotspot script remains executable and reproducible after updates.

---

## Phase G Additions

Phase G introduces a dedicated model-ready benchmark lane for `TotalCrossSection` to avoid fragile generic fixtures.

### Implemented changes

- Added dedicated profiling script:
    - `scripts/dev/profile_total_cross_section_model_ready.jl`
- Script now builds model-ready inputs explicitly:
    - auto-populates `A` via `ensure_quark_params_has_A`
    - computes full effective `K_coeffs` via `calculate_G_from_A` + `calculate_effective_couplings`
- Benchmarks three TCS API surfaces:
    - `total_cross_section`
    - `scan_s_dependence`
    - `calculate_all_total_cross_sections`

### Validation summary

- Added smoke test: `tests/unit/relaxtime/test_total_cross_section_model_ready_fixture_smoke.jl`
- Local baseline (ms/call):
    - `total_cross_section(:uu_to_uu, s0)` = `0.0168`
    - `scan_s_dependence(4 points)` = `0.0529`
    - `calculate_all_total_cross_sections` = `1.0785`
- New smoke test passed (`6/6`).

---

## Usage Examples

### Basic Usage

```julia
using Main.ParameterTypes: QuarkParams, ThermoParams

# Create parameters
q = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
t = ThermoParams(0.15, 0.5, 0.5, 0.0)

# Use in calculations
result = relaxation_times(q, t, K_coeffs; densities=densities)
```

### Conversion Between Formats

```julia
# Struct to NamedTuple
q = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
q_nt = as_namedtuple(q)

# NamedTuple to Struct
q_nt = (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
q = QuarkParams(q_nt)

# Round-trip conversion
q_original = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
q_nt = as_namedtuple(q_original)
q_reconstructed = QuarkParams(q_nt)
# q_original == q_reconstructed
```

### Dynamic Parameter Construction

```julia
# Build parameters programmatically
function create_params(T, μ_B)
    q = QuarkParams(
        m=(u=1.52, d=1.52, s=3.04),
        μ=(u=μ_B/3, d=μ_B/3, s=μ_B/3)
    )
    
    t = ThermoParams(T, 0.5, 0.5, 0.0)
    
    return q, t
end

q, t = create_params(0.15, 0.8)
```

### Parameter Modification

```julia
# Modify parameters using NamedTuple operations
q = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))

# Convert to NamedTuple for modification
q_nt = as_namedtuple(q)

# Modify chemical potentials
q_modified_nt = merge(q_nt, (μ=(u=0.5, d=0.5, s=0.5),))

# Convert back to struct
q_modified = QuarkParams(q_modified_nt)
```

### Error Handling

```julia
# Structs catch errors at construction time
try
    q = QuarkParams((m=(u=1.52, d=1.52, s=3.04),))  # Missing μ
catch e
    println("Error: ", e)
    # Provides clear error message about missing field
end

# NamedTuples catch errors at usage time
q_nt = (m=(u=1.52, d=1.52, s=3.04),)  # Missing μ
try
    result = relaxation_times(q_nt, t, K_coeffs; densities=densities)
catch e
    println("Error: ", e)
    # Error occurs when trying to access μ field
end
```

---

## Type Hierarchy

```
Any
├── QuarkParams
│   ├── m::NamedTuple{(:u, :d, :s), Tuple{Float64, Float64, Float64}}
│   └── μ::NamedTuple{(:u, :d, :s), Tuple{Float64, Float64, Float64}}
└── ThermoParams
    ├── T::Float64
    ├── Φ::Float64
    ├── Φbar::Float64
    └── ξ::Float64
```

---

## Performance Considerations

### Zero Overhead Abstraction

The dual interface pattern has **zero runtime overhead** because:

1. **Inlined normalization**: `@inline` ensures helpers are compiled away
2. **Type specialization**: Julia compiles separate versions for struct and NamedTuple inputs
3. **Constant propagation**: Type checks are resolved at compile time

### Benchmarks

```julia
using BenchmarkTools

q_struct = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
q_nt = as_namedtuple(q_struct)

# Both versions have identical performance
@btime relaxation_times($q_struct, $t_struct, $K; densities=$densities)
@btime relaxation_times($q_nt, $t_nt, $K; densities=$densities)
```

Results show no measurable performance difference between struct and NamedTuple usage.

---

## Phase B Additions

### WorkflowParamAdapters

Phase B introduces a shared workflow adapter module:

- `src/pnjl/workflows/WorkflowParamAdapters.jl`

Provided helpers:

```julia
normalize_quark_params(q)
normalize_thermo_params(t)
as_legacy_inputs(q, t)
```

Compatibility note:
- `as_legacy_inputs(q, t)` is deprecated and kept only for compatibility.
- Prefer `as_relaxtime_inputs(q, t)` in new code.
- The deprecated helper is no longer part of the default recommended surface and should only be used by explicit compatibility paths.

Contract:
- struct input (`QuarkParams` / `ThermoParams`) is converted to legacy-compatible NamedTuple;
- NamedTuple input is passed through unchanged;
- `as_legacy_inputs(q, t)` returns a pair for old interfaces only:
    - `quark_params`
    - `thermo_params`

Used by:
- `src/pnjl/workflows/TransportWorkflow.jl`
- `src/pnjl/workflows/MesonMassWorkflow.jl`

### ScanConfig (structured scan configuration)

Phase B adds scan configuration objects:

- `src/pnjl/scans/ScanConfig.jl`
- `TmuScanConfig`
- `TrhoScanConfig`
- `scan_kwargs(cfg)`

Compatibility pattern:

```julia
run_tmu_scan(; kwargs...)                        # existing API
run_tmu_scan(config::TmuScanConfig; kwargs...)   # new API

run_trho_scan(; kwargs...)                       # existing API
run_trho_scan(config::TrhoScanConfig; kwargs...) # new API
```

Behavior:
- existing kwargs workflows continue to work unchanged;
- config-object style is now supported;
- explicit kwargs override config fields.

---

## Best Practices

### When to Use Structs

✅ **Use structs for:**
- New code and new features
- Public APIs and library interfaces
- Code that benefits from type safety
- Long-lived parameter objects
- Code that needs IDE autocomplete

### When to Use NamedTuples

✅ **Use NamedTuples for:**
- Existing code (no need to migrate)
- Quick prototyping and experimentation
- Dynamic parameter construction
- Interfacing with legacy code
- Inline parameter literals

### Migration Strategy

1. **Don't rush**: Existing NamedTuple code works fine
2. **Migrate incrementally**: Update code as you touch it
3. **Start with new code**: Use structs for all new features
4. **Test thoroughly**: Verify numerical equivalence after migration
5. **Document changes**: Update comments and docstrings

---

## See Also

- [Parameter Struct Migration Guide](../guides/PARAMETER_STRUCT_MIGRATION.md)
- [Usage Examples](../guides/examples/struct_usage_examples.jl)
- [Quick Start Guide](../guides/examples/quick_start_structs.jl)
- [Migration Example](../guides/examples/migration_example.jl)
