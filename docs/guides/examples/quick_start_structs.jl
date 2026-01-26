# Quick Start: Using Parameter Structs
#
# This file shows the essential patterns for using QuarkParams and ThermoParams
# structs in the PNJL model codebase.

using Main.ParameterTypes: QuarkParams, ThermoParams

# =============================================================================
# Pattern 1: Struct Usage (Recommended for New Code)
# =============================================================================

# Create parameters using structs
quark_params = QuarkParams(
    m=(u=1.52, d=1.52, s=3.04),  # Quark masses in fm⁻¹
    μ=(u=0.3, d=0.3, s=0.3)      # Chemical potentials in fm⁻¹
)

thermo_params = ThermoParams(
    0.15,  # Temperature T in fm⁻¹
    0.5,   # Polyakov loop Φ
    0.5,   # Conjugate Polyakov loop Φbar
    0.0    # Anisotropy parameter ξ
)

# Use in calculations
# result = relaxation_times(quark_params, thermo_params, K_coeffs; densities=densities)

# =============================================================================
# Pattern 2: NamedTuple Usage (Backward Compatible)
# =============================================================================

# Create parameters using NamedTuples (old style, still works!)
quark_params_nt = (
    m=(u=1.52, d=1.52, s=3.04),
    μ=(u=0.3, d=0.3, s=0.3)
)

thermo_params_nt = (
    T=0.15,
    Φ=0.5,
    Φbar=0.5,
    ξ=0.0
)

# Use in calculations (produces identical results)
# result = relaxation_times(quark_params_nt, thermo_params_nt, K_coeffs; densities=densities)

# =============================================================================
# Pattern 3: Mixed Usage
# =============================================================================

# You can mix structs and NamedTuples in the same call
# result = relaxation_times(
#     quark_params,      # Struct
#     thermo_params_nt,  # NamedTuple
#     K_coeffs;
#     densities=densities
# )

# =============================================================================
# Pattern 4: Converting Between Formats
# =============================================================================

using Main.ParameterTypes: as_namedtuple

# Convert struct to NamedTuple
quark_nt = as_namedtuple(quark_params)

# Convert NamedTuple to struct
quark_struct = QuarkParams(quark_nt)

# =============================================================================
# Why Use Structs?
# =============================================================================

# ✅ Type safety: Catch missing fields at construction time
# ✅ Better IDE support: Autocomplete and type hints
# ✅ Clearer intent: Explicit parameter types
# ✅ Validation: Constructors can validate field values
# ✅ Documentation: Self-documenting parameter structure

# =============================================================================
# When to Use NamedTuples?
# =============================================================================

# ✅ Existing code: No need to modify working scripts
# ✅ Quick prototyping: Faster to write inline
# ✅ Dynamic construction: When building parameters programmatically
# ✅ Legacy interfaces: When interfacing with older code

# =============================================================================
# Key Takeaway
# =============================================================================

# Both formats work identically - use structs for new code, but don't worry
# about updating existing code unless you want the benefits of type safety.
