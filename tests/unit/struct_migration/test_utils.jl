"""
Test utilities for parameter struct migration tests.

Provides:
- Random parameter generators for property-based testing
- Helper functions for comparing results
- Common test fixtures
"""

using Test
using Supposition

# Load ParameterTypes module
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, joinpath(@__DIR__, "../../../src/ParameterTypes.jl"))
end

using Main.ParameterTypes: QuarkParams, ThermoParams, as_namedtuple

# ============================================================================
# Random Parameter Generators for Property-Based Testing
# ============================================================================

# Note: Supposition.jl uses inline generators in @check functions
# These helper functions provide common test fixtures and comparison utilities

# ============================================================================
# Comparison Helpers
# ============================================================================

"""
    approx_equal(x, y; rtol=1e-12, atol=1e-14)

Compare two values with appropriate tolerance for floating-point comparison.
Uses relative tolerance for large values and absolute tolerance for small values.
"""
function approx_equal(x::Real, y::Real; rtol=1e-12, atol=1e-14)
    return isapprox(x, y, rtol=rtol, atol=atol)
end

function approx_equal(x::NamedTuple, y::NamedTuple; rtol=1e-12, atol=1e-14)
    keys(x) == keys(y) || return false
    return all(approx_equal(getfield(x, k), getfield(y, k); rtol=rtol, atol=atol) for k in keys(x))
end

function approx_equal(x::Tuple, y::Tuple; rtol=1e-12, atol=1e-14)
    length(x) == length(y) || return false
    return all(approx_equal(x[i], y[i]; rtol=rtol, atol=atol) for i in 1:length(x))
end

# ============================================================================
# Common Test Fixtures
# ============================================================================

"""
    get_test_quark_params()

Return standard test QuarkParams for unit tests.
"""
function get_test_quark_params()
    return (
        m = (u=1.52, d=1.52, s=3.04),
        μ = (u=0.3, d=0.3, s=0.3)
    )
end

"""
    get_test_thermo_params()

Return standard test ThermoParams for unit tests.
"""
function get_test_thermo_params()
    return (
        T = 0.15,
        Φ = 0.5,
        Φbar = 0.5,
        ξ = 0.0
    )
end

"""
    get_test_K_coeffs()

Return standard test coupling coefficients for unit tests.
This is a placeholder - actual values should come from the codebase.
"""
function get_test_K_coeffs()
    # These are placeholder values - actual K_coeffs structure
    # should be determined from the existing codebase
    return (
        K_pi = 1.0,
        K_K = 1.0,
        K_eta = 1.0,
        K_sigma = 1.0,
        K_kappa = 1.0,
        K_a0 = 1.0
    )
end

# ============================================================================
# Exports
# ============================================================================

export QuarkParams, ThermoParams, as_namedtuple
export approx_equal
export get_test_quark_params, get_test_thermo_params, get_test_K_coeffs

