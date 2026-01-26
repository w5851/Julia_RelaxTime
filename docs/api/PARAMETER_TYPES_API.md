# Parameter Types API Documentation

This document provides comprehensive API documentation for the `QuarkParams` and `ThermoParams` structs and their associated helper functions.

## Table of Contents

- [QuarkParams](#quarkparams)
- [ThermoParams](#thermoparams)
- [Helper Functions](#helper-functions)
- [Normalization Helpers](#normalization-helpers)
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
