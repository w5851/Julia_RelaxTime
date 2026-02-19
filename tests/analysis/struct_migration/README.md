# Parameter Struct Migration Tests

This directory contains tests for the parameter struct migration project, which migrates the PNJL model codebase from mixed NamedTuple/struct parameter passing to a unified struct-based approach.

## Directory Structure

- `test_utils.jl` - Shared test utilities for generating random valid parameters
- `test_parameter_types.jl` - Tests for QuarkParams and ThermoParams struct construction and conversion
- `test_particle_symbols_struct.jl` - Tests for ParticleSymbols module with struct support
- `test_differential_cross_section_struct.jl` - Tests for DifferentialCrossSection module with struct support
- `test_scattering_amplitude_struct.jl` - Tests for ScatteringAmplitude module with struct support
- `test_total_propagator_struct.jl` - Tests for TotalPropagator module with struct support
- `test_total_cross_section_struct.jl` - Tests for TotalCrossSection module with struct support
- `test_average_scattering_rate_struct.jl` - Tests for AverageScatteringRate module with struct support
- `test_relaxation_time_struct.jl` - Tests for RelaxationTime module with struct support
- `test_integration.jl` - Integration tests for full call chain with struct parameters
- `baseline_results.txt` - Baseline test results from existing test suite for comparison

## Test Strategy

The migration uses a dual testing approach:

### Unit Tests
- Verify specific examples with known correct outputs
- Test edge cases (zero masses, zero temperature, etc.)
- Test error conditions (missing fields, invalid inputs)
- Verify backward compatibility with existing test cases

### Property-Based Tests
- Verify struct-NamedTuple equivalence across random inputs
- Test field extraction correctness with random parameters
- Verify conversion round-trip preservation
- Test numerical stability across parameter ranges

All property-based tests use Supposition.jl with a minimum of 100 iterations per test.

## Running Tests

Run all struct migration tests:
```julia
julia --project=. -e 'using Test; include("tests/unit/struct_migration/runtests.jl")'
```

Run a specific test file:
```julia
julia --project=. -e 'using Test; include("tests/unit/struct_migration/test_parameter_types.jl")'
```

## Requirements Coverage

These tests validate requirements 1-12 from the parameter struct migration specification:
- Requirements 1-7: Dual interface support for all modules
- Requirement 8: Backward compatibility preservation
- Requirement 9: Internal consistency and validation
- Requirement 10: Performance preservation
- Requirement 11: Test coverage for struct equivalence
- Requirement 12: Documentation and usage guidance
