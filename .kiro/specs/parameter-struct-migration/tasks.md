# Implementation Plan: Parameter Struct Migration

## Overview

This implementation plan breaks down the parameter struct migration into discrete, incremental tasks. The migration follows a bottom-up approach, starting with utility modules and working up to the main RelaxationTime module. Each task builds on previous work and includes testing to catch errors early.

The plan prioritizes:
1. **Incremental progress**: Each task produces working, tested code
2. **Early validation**: Tests verify correctness at each step
3. **Minimal disruption**: Existing code continues to work throughout
4. **Clear dependencies**: Tasks build logically on each other

## Tasks

- [x] 1. Set up testing infrastructure and baseline
  - Create test directory structure for struct migration tests
  - Run existing test suite and record baseline results for comparison
  - Set up Supposition.jl for property-based testing
  - Create test utilities for generating random valid parameters
  - _Requirements: 8.1, 11.1-11.7_

- [ ] 2. Update ParticleSymbols utility module
  - [x] 2.1 Add struct support to ParticleSymbols.jl
    - Import QuarkParams from Main.ParameterTypes
    - Add normalization helper `_nt_quark(q)`
    - Update `get_mass` signature to accept `Union{NamedTuple, QuarkParams}`
    - Update `get_chemical_potential` signature to accept `Union{NamedTuple, QuarkParams}`
    - Update `get_quark_masses_for_process` signature to accept `Union{NamedTuple, QuarkParams}`
    - Add normalization calls at function entry points
    - Update docstrings to document both parameter formats
    - _Requirements: 7.1, 7.2, 7.3, 7.5_
  
  - [x] 2.2 Write property tests for ParticleSymbols
    - **Property 2: Field Extraction Correctness**
    - **Validates: Requirements 7.1, 7.2, 7.3**
    - Test `get_mass` with random QuarkParams structs vs NamedTuples
    - Test `get_chemical_potential` with random QuarkParams structs vs NamedTuples
    - Test `get_quark_masses_for_process` with random parameters
    - Verify extracted values match struct fields exactly
  
  - [x] 2.3 Write unit tests for ParticleSymbols edge cases
    - Test error handling for invalid particle symbols
    - Test all particle types (:u, :d, :s, :ubar, :dbar, :sbar)
    - Test all scattering processes
    - _Requirements: 7.1, 7.2, 7.3_

- [ ] 3. Update DifferentialCrossSection module
  - [x] 3.1 Add struct support to DifferentialCrossSection.jl
    - Import QuarkParams and ThermoParams from Main.ParameterTypes
    - Add normalization helpers `_nt_quark` and `_nt_thermo`
    - Update `differential_cross_section` signature to accept Union types
    - Add normalization calls at function entry
    - Update docstrings
    - _Requirements: 5.1, 5.2, 5.3_
  
  - [x] 3.2 Write property test for differential cross-section equivalence
    - **Property 1: Struct-NamedTuple Equivalence (DifferentialCrossSection)**
    - **Validates: Requirements 5.2, 11.5**
    - Generate random valid parameters (masses, chemical potentials, T, Φ, Φbar)
    - Generate random kinematic variables (s, t within physical bounds)
    - Compute differential cross-section with struct and NamedTuple parameters
    - Verify numerical equivalence within tolerance (rtol=1e-12)
  
  - [x] 3.3 Write unit tests for differential cross-section edge cases
    - Test threshold behavior (s near kinematic threshold)
    - Test high-energy behavior
    - Test different scattering processes
    - _Requirements: 5.1, 5.5_

- [ ] 4. Update ScatteringAmplitude module
  - [x] 4.1 Add struct support to ScatteringAmplitude.jl
    - Import QuarkParams and ThermoParams from Main.ParameterTypes
    - Add normalization helpers
    - Update `scattering_amplitude_squared` signature to accept Union types
    - Add normalization calls at function entry
    - Update docstrings
    - _Requirements: 4.1, 4.2, 4.3_
  
  - [x] 4.2 Write property test for scattering amplitude equivalence
    - **Property 1: Struct-NamedTuple Equivalence (ScatteringAmplitude)**
    - **Validates: Requirements 4.2, 11.4**
    - Generate random valid parameters
    - Generate random kinematic variables (s, t)
    - Compute amplitude squared with struct and NamedTuple parameters
    - Verify numerical equivalence within tolerance
    - Verify amplitude is non-negative (since it's squared)
  
  - [x] 4.3 Write unit tests for amplitude calculations
    - Test different scattering processes
    - Test amplitude with different K_coeffs values
    - Test that amplitude depends correctly on coupling strengths
    - _Requirements: 4.4, 4.5_

- [ ] 5. Update TotalPropagator module
  - [x] 5.1 Add struct support to TotalPropagator.jl
    - Import QuarkParams and ThermoParams from Main.ParameterTypes
    - Add normalization helpers
    - Update `total_propagator_simple` signature to accept Union types
    - Update `total_propagator_mixed` signature to accept Union types
    - Update `total_propagator_auto` signature to accept Union types
    - Update `calculate_cms_momentum` signature to accept Union types
    - Update `calculate_all_propagators` signature to accept Union types
    - Add normalization calls at all function entry points
    - Update docstrings
    - _Requirements: 6.1, 6.2, 6.3, 6.4, 6.5_
  
  - [x] 5.2 Write property tests for propagator equivalence
    - **Property 1: Struct-NamedTuple Equivalence (TotalPropagator)**
    - **Validates: Requirements 6.2, 11.6**
    - Test `total_propagator_simple` with random parameters
    - Test `total_propagator_mixed` with random parameters
    - Test `calculate_cms_momentum` with random kinematic variables
    - Verify numerical equivalence for all functions
  
  - [x] 5.3 Write unit tests for propagator calculations
    - Test different scattering processes and channels
    - Test different meson combinations
    - Test kinematic edge cases
    - _Requirements: 6.1, 6.3, 6.4_

- [x] 6. Update TotalCrossSection module
  - [x] 6.1 Add struct support to TotalCrossSection.jl
    - Import QuarkParams and ThermoParams from Main.ParameterTypes
    - Add normalization helpers
    - Update `total_cross_section` signature to accept Union types
    - Update `calculate_all_total_cross_sections` signature to accept Union types
    - Update `scan_s_dependence` signature to accept Union types
    - Update helper functions (`calculate_t_bounds`, `calculate_final_state_energies`, etc.)
    - Add normalization calls at function entry points
    - Update docstrings
    - _Requirements: 3.1, 3.2, 3.3, 3.4, 3.5_
    - _Status: COMPLETED - Module already has full struct support_
  
  - [x] 6.2 Write property test for total cross-section equivalence
    - **Property 1: Struct-NamedTuple Equivalence (TotalCrossSection)**
    - **Validates: Requirements 3.2, 11.3**
    - Generate random valid parameters
    - Generate random s values above threshold
    - Compute total cross-section with struct and NamedTuple parameters
    - Verify numerical equivalence within tolerance
    - Verify cross-section is non-negative
  
  - [x] 6.3 Write unit tests for cross-section calculations
    - Test `calculate_all_total_cross_sections` with struct parameters
    - Test `scan_s_dependence` produces smooth curves
    - Test threshold behavior
    - Test different scattering processes
    - _Requirements: 3.4, 3.5_

- [x] 7. Update AverageScatteringRate module
  - [x] 7.1 Add struct support to AverageScatteringRate.jl
    - Import QuarkParams and ThermoParams from Main.ParameterTypes
    - Add normalization helpers `_nt_quark` and `_nt_thermo`
    - Update `average_scattering_rate` signature to accept Union types
    - Update `build_w0cdf_pchip_cache` signature to accept Union types
    - Update `precompute_cross_section!` signature to accept Union types
    - Update helper functions (`get_mass`, `get_mu`) to accept Union types
    - Add normalization calls at function entry points
    - Update docstrings
    - _Requirements: 2.1, 2.2, 2.3, 2.4, 2.5_
  
  - [x] 7.2 Write property test for average scattering rate equivalence
    - **Property 1: Struct-NamedTuple Equivalence (AverageScatteringRate)**
    - **Validates: Requirements 2.2, 11.2**
    - Generate random valid parameters
    - Compute average scattering rate with struct and NamedTuple parameters
    - Verify numerical equivalence within tolerance
    - Verify rate is non-negative
  
  - [x] 7.3 Write unit tests for scattering rate calculations
    - Test `build_w0cdf_pchip_cache` with struct parameters
    - Test cache construction and interpolation
    - Test different scattering processes
    - Test helper functions `get_mass` and `get_mu` with structs
    - _Requirements: 2.4, 2.5_

- [x] 8. Update RelaxationTime main module
  - [x] 8.1 Add struct support to RelaxationTime.jl
    - Verify QuarkParams and ThermoParams are imported (already done)
    - Verify normalization helpers exist (already done)
    - Update `relaxation_times` signature to accept Union types
    - Update `compute_average_rates` signature to accept Union types
    - Add normalization calls at function entry points
    - Verify `ensure_quark_params_has_A` works with structs (already does)
    - Update docstrings
    - _Requirements: 1.1, 1.2, 1.3, 1.4_
  
  - [x] 8.2 Write property test for relaxation times equivalence
    - **Property 1: Struct-NamedTuple Equivalence (RelaxationTime)**
    - **Validates: Requirements 1.2, 11.1**
    - Generate random valid parameters
    - Generate random valid densities
    - Compute relaxation times with struct and NamedTuple parameters
    - Verify numerical equivalence for tau, tau_inv, and rates
    - Verify all relaxation times are positive and finite
  
  - [x] 8.3 Write unit tests for relaxation time calculations
    - Test `compute_average_rates` with struct parameters
    - Test `relaxation_rates` with struct-derived densities
    - Test different scattering process combinations
    - Test with pre-computed rates vs computing from scratch
    - _Requirements: 1.3, 1.4_

- [x] 9. Integration testing and validation
  - [x] 9.1 Write integration tests for full call chain
    - Test that struct parameters flow correctly through entire chain
    - Test RelaxationTime → AverageScatteringRate → TotalCrossSection → ScatteringAmplitude
    - Test with different parameter combinations
    - Verify all intermediate results are consistent
    - _Requirements: 5.4, 8.3_
  
  - [x] 9.2 Write property test for conversion round-trip
    - **Property 3: Conversion Round-Trip Preservation**
    - **Validates: Requirements 9.5**
    - Generate random QuarkParams and ThermoParams structs
    - Convert to NamedTuple using `as_namedtuple`
    - Reconstruct structs from NamedTuples
    - Verify all field values are preserved exactly
  
  - [x] 9.3 Write property test for ensure_quark_params_has_A
    - **Property 5: Extended QuarkParams with A Field**
    - **Validates: Requirements 9.4**
    - Test with QuarkParams that has A field (should preserve it)
    - Test with QuarkParams without A field (should compute it)
    - Test with both struct and NamedTuple inputs
    - Verify output structure is correct

- [x] 10. Backward compatibility verification
  - [x] 10.1 Run existing test suite without modifications
    - **Property 4: Backward Compatibility Preservation**
    - **Validates: Requirements 8.1, 8.2**
    - Run all existing tests in `tests/unit/relaxtime/`
    - Run all existing tests in `tests/unit/pnjl/`
    - Verify all tests that passed before still pass (309 tests)
    - Verify the 1 pre-existing error remains unchanged
    - Verify no numerical results have changed
    - Document any test failures for investigation
  
  - [x] 10.2 Test mixed struct/NamedTuple usage patterns
    - Test calling functions with struct quark_params and NamedTuple thermo_params
    - Test calling functions with NamedTuple quark_params and struct thermo_params
    - Test workflows that mix both representations
    - Verify all combinations work correctly
    - _Requirements: 8.3_

- [x] 11. Documentation and examples
  - [x] 11.1 Update module-level documentation
    - Add explanation of dual interface pattern to each module
    - Document the internal normalization strategy
    - Explain when to use structs vs NamedTuples
    - Add migration guide for existing code
    - _Requirements: 12.1, 12.2, 12.3, 12.4_
  
  - [x] 11.2 Create usage examples
    - Create example showing struct usage (recommended pattern)
    - Create example showing NamedTuple usage (backward compatible)
    - Create example showing mixed usage
    - Create example of converting existing code to use structs
    - _Requirements: 12.3, 12.4_
  
  - [x] 11.3 Update API documentation
    - Update function docstrings to show both parameter formats
    - Add examples to docstrings demonstrating struct usage
    - Document helper functions (`_nt_quark`, `_nt_thermo`)
    - Update README with migration information
    - _Requirements: 12.1, 12.5_

- [x] 12. Performance validation
  - [x] 12.1 Create performance benchmarks
    - Benchmark `relaxation_times` with struct vs NamedTuple parameters
    - Benchmark `average_scattering_rate` with struct vs NamedTuple parameters
    - Benchmark `total_cross_section` with struct vs NamedTuple parameters
    - Benchmark full call chain with both parameter types
    - _Requirements: 10.1, 10.2, 10.3_
  
  - [x] 12.2 Analyze and document performance results
    - Compare struct vs NamedTuple performance
    - Verify performance is within 5% tolerance
    - Document any performance differences
    - Investigate and optimize if needed
    - _Requirements: 10.1, 10.2, 10.3_

- [x] 13. Final validation and cleanup
  - [x] 13.1 Run complete test suite
    - Run all unit tests
    - Run all property-based tests
    - Run all integration tests
    - Run backward compatibility tests
    - Verify all tests pass
  
  - [x] 13.2 Code review and cleanup
    - Review all modified files for consistency
    - Verify all functions have updated docstrings
    - Verify all normalization helpers are inlined
    - Remove any debug code or comments
    - Verify code style consistency
  
  - [x] 13.3 Update migration documentation
    - Update `docs/dev/active/parameter_structs_migration.md`
    - Mark Phase A as complete
    - Document any lessons learned
    - Document any deviations from original plan
    - Prepare notes for Phase B planning

## Notes

- Each task references specific requirements for traceability
- Tasks are ordered to minimize dependencies and enable incremental progress
- Testing tasks are placed immediately after implementation tasks to catch errors early
- The plan assumes Julia language and existing project structure
- Property-based tests use Supposition.jl with minimum 100 iterations
- Unit tests focus on edge cases and specific examples
- Integration tests verify the full call chain works correctly
- Backward compatibility is verified by running existing tests unchanged
