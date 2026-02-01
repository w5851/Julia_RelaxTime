# Implementation Plan: QCD Model Architecture Refactoring

## Overview

This implementation plan breaks down the QCD model architecture refactoring into discrete, actionable coding tasks. The refactoring transforms the monolithic PNJL implementation into a polymorphic, extensible system supporting multiple QCD models (PNJL, rPNJL, and future models) while maintaining backward compatibility and numerical precision.

The implementation follows a phased approach:
1. **Foundation**: Core interfaces, exceptions, and base utilities
2. **PNJL Refactoring**: Extract existing PNJL into new architecture with regression testing
3. **rPNJL Implementation**: Add new rPNJL model with eight-quark terms and Vandermonde potential
4. **Integration**: Factory pattern, compatibility layer, and final validation

## Tasks

- [ ] 1. Set up core infrastructure and interfaces
  - [ ] 1.1 Create core interface definitions
    - Create `src/core/interface.jl` with AbstractQCDModel and domain-specific abstract types (AbstractIsotropicModel, AbstractAnisotropicModel, AbstractMagneticModel)
    - Define required method signatures: vacuum_contribution, thermal_contribution, polyakov_potential, dispersion_relation, calculate_mass_vec, calculate_chiral
    - Add descriptive error messages for unimplemented methods
    - _Requirements: 1.1, 1.2, 1.3, 1.4, 1.5_
  
  - [ ] 1.2 Create exception hierarchy
    - Create `src/core/exceptions.jl` with custom exception types: QCDModelException (base), ModelParameterError, NumericalConvergenceError, PhysicalConstraintError, ModelConfigurationError, NumericalInstabilityError
    - Implement Base.showerror for each exception type with descriptive messages including problem description, context, and suggested actions
    - _Requirements: 13.1, 13.2_
  
  - [ ]* 1.3 Write unit tests for exception hierarchy
    - Test that each exception type can be constructed and thrown
    - Test that error messages include all required information (parameter, value, constraint, suggestion)
    - Test exception inheritance hierarchy
    - _Requirements: 13.7_

- [ ] 2. Implement base utilities layer
  - [ ] 2.1 Create numerical integration module
    - Create `src/models/base/integrals.jl` with integrate_momentum and vacuum_integral functions
    - Implement Gauss-Legendre quadrature with configurable node count
    - Add fallback strategy: if integration fails with n_nodes, retry with 2*n_nodes
    - Accept energy function callbacks for flexibility
    - _Requirements: 2.1, 13.3_
  
  - [ ] 2.2 Create Polyakov loop calculation module
    - Create `src/models/base/polyakov.jl` with polyakov_log_combination, standard_polyakov_potential, and vandermonde_polyakov_potential functions
    - Implement temperature-dependent b2(T) calculation for Vandermonde terms
    - Add input validation for Polyakov loop range [0, 1]
    - _Requirements: 2.2, 4.4, 4.5_
  
  - [ ] 2.3 Create safe math operations module
    - Create `src/models/base/safe_math.jl` with safe_log, safe_vandermonde, and safe_divide functions
    - Implement safe_log with fallback for zero/negative inputs
    - Implement safe_vandermonde with warning logging for negative discriminants but continue calculation
    - Implement safe_divide with zero-denominator checking
    - _Requirements: 2.3, 5.1, 5.2, 5.3, 5.4, 4.6_
  
  - [ ] 2.4 Create distribution functions module
    - Create `src/models/base/distributions.jl` with quark_distribution and antiquark_distribution functions
    - Implement Fermi-Dirac distributions with Polyakov loop modifications
    - Handle zero-temperature limit correctly
    - _Requirements: 2.4_
  
  - [ ] 2.5 Create parameter validation module
    - Create `src/models/base/validation.jl` with validate_parameters and validate_physical_constraint functions
    - Define STANDARD_CONSTRAINTS dictionary with ranges for T, mu, phi, M
    - Implement constraint checking that throws ModelParameterError with descriptive messages
    - _Requirements: 14.1, 14.2, 14.3, 14.4_
  
  - [ ] 2.6 Create mass and chiral condensate module
    - Create `src/models/base/mass_chiral.jl` with solve_mass_gap and calculate_chiral_condensate functions
    - Implement iterative mass gap solver with convergence checking
    - Add fallback strategies: relax tolerance, increase max iterations
    - Throw NumericalConvergenceError if all strategies fail
    - _Requirements: 2.5, 13.4_
  
  - [ ]* 2.7 Write unit tests for base utilities
    - Test integrate_momentum with known analytical integrals
    - Test polyakov_log_combination with edge cases (phi=0, phi=1, phi=0.5)
    - Test safe_log with zero, negative, and positive inputs
    - Test safe_vandermonde with negative discriminants
    - Test distribution functions return values in [0, 1]
    - Test parameter validation with valid and invalid inputs
    - Test mass gap solver convergence and fallback strategies
    - _Requirements: 8.1_
  
  - [ ]* 2.8 Write property test for numerical integration (Property 2)
    - **Property 2: Numerical integration produces finite results**
    - Generate random valid energy functions, temperatures, and chemical potentials
    - Verify integrate_momentum returns finite Float64 (not NaN or Inf)
    - Run with 100+ iterations
    - **Validates: Requirements 2.1**
  
  - [ ]* 2.9 Write property test for Polyakov log calculations (Property 3)
    - **Property 3: Polyakov log calculations respect physical constraints**
    - Generate random phi and phi_bar values in [0, 1]
    - Verify polyakov_log_combination returns finite values
    - Test edge cases: phi=0, phi=1, phi_bar=0, phi_bar=1
    - Run with 100+ iterations
    - **Validates: Requirements 2.2**
  
  - [ ]* 2.10 Write property test for safe math functions (Property 4)
    - **Property 4: Safe math functions handle edge cases without crashing**
    - Generate random inputs including edge cases (zero, negative, very large, very small)
    - Verify safe_log, safe_vandermonde, safe_divide either return valid results or throw descriptive exceptions
    - Verify no unexpected NaN or crashes
    - Run with 100+ iterations
    - **Validates: Requirements 2.3, 5.1, 5.2**
  
  - [ ]* 2.11 Write property test for distribution functions (Property 5)
    - **Property 5: Distribution functions return values in valid range**
    - Generate random valid energy, temperature, chemical potential, and Polyakov loop values
    - Verify quark_distribution and antiquark_distribution return values in [0, 1]
    - Run with 100+ iterations
    - **Validates: Requirements 2.4**

- [ ] 3. Checkpoint - Base utilities validation
  - Ensure all base utilities unit tests pass
  - Ensure all base utilities property tests pass (Properties 2-5)
  - Ask the user if questions arise

- [ ] 4. Refactor PNJL model into new architecture
  - [ ] 4.1 Create PNJL model type and structure
    - Create `src/models/isotropic/types.jl` with AbstractIsotropicModel definition
    - Create `src/models/isotropic/pnjl.jl` with PNJLModel struct inheriting from AbstractIsotropicModel
    - Define fields: params (QuarkParams), thermo_params (ThermoParams), cutoff
    - Add constructor with parameter validation
    - _Requirements: 3.1, 12.1_
  
  - [ ] 4.2 Implement PNJL interface methods
    - Implement vacuum_contribution using vacuum_integral from base utilities
    - Implement thermal_contribution using integrate_momentum and solve_mass_gap
    - Implement polyakov_potential using standard_polyakov_potential
    - Implement dispersion_relation as sqrt(p^2 + M^2)
    - Implement calculate_mass_vec by solving mass gap for each flavor
    - Implement calculate_chiral using calculate_chiral_condensate
    - _Requirements: 3.3, 3.4, 2.5_
  
  - [ ] 4.3 Create parameter data structures
    - Create QuarkParams struct with fields: n_flavors, current_masses, coupling_G, coupling_K, cutoff
    - Create ThermoParams struct with fields: T_range, μ_range, n_T_points, n_μ_points, Polyakov potential parameters (a0-a3, b3-b4)
    - Add constructors with validation using @assert
    - Implement immutable design for thread safety
    - _Requirements: 12.1, 14.6_
  
  - [ ] 4.4 Implement parameter format conversion
    - Add QuarkParams constructor accepting NamedTuple
    - Add QuarkParams constructor accepting Dict{Symbol, Any}
    - Add Base.Dict method to convert QuarkParams to Dict
    - Implement same conversions for ThermoParams
    - _Requirements: 12.2, 12.3, 12.5_
  
  - [ ] 4.5 Create regression test suite for PNJL
    - Create `tests/regression/pnjl_refactoring/test_thermodynamics.jl`
    - Load baseline data from original PNJL implementation
    - Compare vacuum_contribution, thermal_contribution, polyakov_potential for multiple (T, μ) points
    - Compare constituent masses and chiral condensates
    - Assert relative error < 1e-12 for all quantities
    - _Requirements: 3.2, 8.3_
  
  - [ ]* 4.6 Write unit tests for PNJL model
    - Test PNJLModel construction with valid and invalid parameters
    - Test each interface method returns expected types
    - Test specific thermodynamic calculations with known results
    - Test edge cases: T=0, μ=0, extreme values
    - _Requirements: 8.1_
  
  - [ ]* 4.7 Write property test for PNJL numerical equivalence (Property 6)
    - **Property 6: Refactored PNJL produces numerically equivalent results**
    - Generate random valid temperatures and chemical potentials
    - Calculate thermodynamic quantities with refactored PNJLModel
    - Compare with original implementation results
    - Verify relative error < 1e-12 for all quantities
    - Run with 100+ iterations
    - **Validates: Requirements 3.2**
  
  - [ ]* 4.8 Write property test for interface compliance (Property 1)
    - **Property 1: All concrete models implement required interface methods**
    - For PNJLModel instance, call each required interface method
    - Verify no MethodError is thrown
    - Verify each method returns appropriate type
    - Run with 100+ iterations with random parameters
    - **Validates: Requirements 1.3**

- [ ] 5. Implement backward compatibility layer
  - [ ] 5.1 Create legacy API wrapper
    - Create `src/compatibility/legacy_api.jl` with legacy function wrappers
    - Implement calculate_thermodynamics_legacy that converts old parameters to new format
    - Implement pnjl_thermodynamics wrapper using default parameters
    - Implement convert_legacy_params to map old parameter names to new QuarkParams
    - _Requirements: 6.1, 6.2, 6.3_
  
  - [ ] 5.2 Maintain public exports
    - Ensure all existing public functions are exported
    - Add compatibility layer functions to exports
    - Document which functions are legacy vs new API
    - _Requirements: 6.4_
  
  - [ ]* 5.3 Write unit tests for compatibility layer
    - Test legacy API functions produce same results as new API
    - Test parameter conversion handles old and new formats
    - Test all exported functions are accessible
    - _Requirements: 8.1_
  
  - [ ]* 5.4 Write property test for legacy API equivalence (Property 7)
    - **Property 7: Legacy API produces equivalent results through compatibility layer**
    - Generate random valid legacy API calls with old parameter format
    - Compare results with original system
    - Verify relative error < 1e-12
    - Run with 100+ iterations
    - **Validates: Requirements 3.5, 6.2, 6.5**

- [ ] 6. Checkpoint - PNJL refactoring validation
  - Ensure all PNJL unit tests pass
  - Ensure all regression tests pass with <1e-12 error
  - Ensure Properties 1, 6, 7 pass
  - Verify backward compatibility with existing scripts
  - Ask the user if questions arise

- [ ] 7. Implement rPNJL model
  - [ ] 7.1 Create rPNJL model type and structure
    - Create `src/models/anisotropic/types.jl` with AbstractAnisotropicModel definition
    - Create `src/models/anisotropic/rpnjl.jl` with RPNJLModel struct inheriting from AbstractAnisotropicModel
    - Define fields: params, thermo_params, cutoff, g1, g2, kappa, T0
    - Add constructor with validation for rPNJL-specific parameters (g1, g2, kappa, T0 ranges)
    - _Requirements: 4.1, 4.2_
  
  - [ ] 7.2 Implement rPNJL interface methods
    - Implement vacuum_contribution (same as PNJL)
    - Implement thermal_contribution including eight-quark term contribution
    - Implement polyakov_potential using vandermonde_polyakov_potential
    - Implement dispersion_relation (same as PNJL)
    - Implement calculate_mass_vec using solve_mass_gap_rpnjl with eight-quark corrections
    - Implement calculate_chiral (same as PNJL)
    - _Requirements: 4.3, 4.4_
  
  - [ ] 7.3 Implement eight-quark interaction terms
    - Implement calculate_eight_quark_contribution function (equation 3.31)
    - Calculate eight-quark term: -g1 * σ^4 - g2 * σ^2 * Δ^2
    - Implement calculate_eight_quark_mass_correction for mass gap equation
    - _Requirements: 4.3_
  
  - [ ] 7.4 Implement Vandermonde-modified Polyakov potential
    - Implement calculate_b2_temperature function for b2(T) temperature dependence
    - Implement Vandermonde term calculation (equations 3.27-3.29)
    - Use safe_vandermonde for negative discriminant handling with warning logging
    - Combine with standard Polyakov potential
    - _Requirements: 4.4, 4.5, 4.6_
  
  - [ ] 7.5 Implement rPNJL mass gap solver
    - Implement solve_mass_gap_rpnjl with eight-quark corrections
    - Combine standard PNJL mass term with eight-quark correction
    - Add convergence checking and fallback strategies
    - Throw NumericalConvergenceError if solver fails
    - _Requirements: 4.3, 13.4_
  
  - [ ]* 7.6 Write unit tests for rPNJL model
    - Test RPNJLModel construction with valid and invalid parameters
    - Test eight-quark term calculation with known values
    - Test Vandermonde term calculation with various temperatures
    - Test b2(T) temperature dependence
    - Test mass gap solver convergence
    - Test edge cases and error conditions
    - _Requirements: 8.1, 8.5_
  
  - [ ]* 7.7 Write property test for rPNJL eight-quark terms (Property 8)
    - **Property 8: rPNJL mass calculation includes eight-quark terms**
    - Generate random valid temperatures, chemical potentials, and rPNJL parameters (g1 > 0 or g2 > 0)
    - Calculate constituent mass with RPNJLModel
    - Calculate constituent mass with PNJLModel (same base parameters)
    - Verify masses differ in a way consistent with repulsive interaction (rPNJL mass reflects eight-quark contribution)
    - Run with 100+ iterations
    - **Validates: Requirements 4.3**
  
  - [ ]* 7.8 Write property test for rPNJL Vandermonde terms (Property 9)
    - **Property 9: rPNJL Polyakov potential includes Vandermonde terms**
    - Generate random valid temperatures and Vandermonde parameters (kappa, T0)
    - Calculate Polyakov potential with RPNJLModel
    - Calculate standard polynomial Polyakov potential
    - Verify potentials differ in a way consistent with Vandermonde contribution
    - Run with 100+ iterations
    - **Validates: Requirements 4.4**
  
  - [ ]* 7.9 Write property test for Vandermonde temperature dependence (Property 10)
    - **Property 10: Vandermonde term exhibits temperature dependence**
    - Generate random pairs of different temperatures (T1 ≠ T2) and fixed Polyakov loop value
    - Calculate Vandermonde term at both temperatures
    - Verify terms produce different values (not temperature-independent)
    - Run with 100+ iterations
    - **Validates: Requirements 4.5**

- [ ] 8. Implement factory pattern and model creation
  - [ ] 8.1 Create model factory
    - Create `src/models/factory.jl` with create_model function
    - Accept model_name symbol (:PNJL, :rPNJL, :MagneticPNJL)
    - Dispatch to model-specific creation functions (create_pnjl_model, create_rpnjl_model)
    - Throw ModelConfigurationError for unknown model names with list of available models
    - _Requirements: 9.1, 9.2, 9.3, 9.5_
  
  - [ ] 8.2 Implement model-specific factory functions
    - Implement create_pnjl_model accepting params, thermo_params, cutoff
    - Implement create_rpnjl_model accepting params, thermo_params, cutoff, g1, g2, kappa, T0
    - Add keyword argument support for flexible parameter passing
    - _Requirements: 9.4_
  
  - [ ]* 8.3 Write unit tests for factory pattern
    - Test create_model with :PNJL returns PNJLModel instance
    - Test create_model with :rPNJL returns RPNJLModel instance
    - Test create_model with unknown name throws ModelConfigurationError
    - Test factory accepts keyword arguments correctly
    - _Requirements: 8.1_

- [ ] 9. Implement parameter management and validation
  - [ ]* 9.1 Write property test for parameter conversion (Property 11)
    - **Property 11: Parameter format conversion preserves values**
    - Generate random valid parameter sets
    - Convert QuarkParams → Dict → QuarkParams
    - Convert QuarkParams → NamedTuple → QuarkParams (if supported)
    - Verify all parameter values preserved (round-trip property)
    - Run with 100+ iterations
    - **Validates: Requirements 12.2, 12.3, 12.5**
  
  - [ ]* 9.2 Write property test for parameter format acceptance (Property 12)
    - **Property 12: Parameter conversion accepts all supported formats**
    - Generate random valid parameter data in each supported format (struct, NamedTuple, Dict{Symbol, Any})
    - Verify conversion utilities successfully convert to target format without exceptions
    - Run with 100+ iterations
    - **Validates: Requirements 6.3**
  
  - [ ]* 9.3 Write property test for parameter validation (Property 13)
    - **Property 13: Invalid parameters trigger validation errors**
    - Generate random invalid parameter values (negative temperature, Polyakov loop outside [0,1], negative mass, etc.)
    - Attempt model construction or calculation
    - Verify ModelParameterError is thrown with descriptive message indicating violated constraint
    - Run with 100+ iterations
    - **Validates: Requirements 14.2, 14.3**
  
  - [ ]* 9.4 Write property test for physical constraints (Property 14)
    - **Property 14: Physical constraints are enforced during calculations**
    - Generate random calculations that may produce constraint violations
    - Verify system validates constraints (Polyakov loop in [0,1], non-negative masses, finite energies)
    - Verify PhysicalConstraintError thrown when violated
    - Run with 100+ iterations
    - **Validates: Requirements 14.4**
  
  - [ ]* 9.5 Write property test for exception messages (Property 15)
    - **Property 15: All exceptions include descriptive messages**
    - Trigger each custom exception type (ModelParameterError, NumericalConvergenceError, PhysicalConstraintError, ModelConfigurationError, NumericalInstabilityError)
    - Verify each exception contains non-empty descriptive message
    - Verify message includes context (what went wrong, parameter values, suggested remediation)
    - Run with 100+ iterations
    - **Validates: Requirements 13.2**

- [ ] 10. Implement logging and monitoring
  - [ ] 10.1 Add structured logging throughout codebase
    - Add @debug logs for detailed internal state (mass gap iterations, integration parameters)
    - Add @info logs for normal operation milestones (model creation, calculation completion, scan start)
    - Add @warn logs for recoverable issues (slow convergence, negative Vandermonde discriminant, performance degradation)
    - Add @error logs for serious errors (parameter validation failures, solver divergence)
    - Use structured format with key-value pairs: @info "Operation" key1=value1 key2=value2
    - _Requirements: 15.1, 15.2, 15.3, 15.4, 15.5, 15.6_
  
  - [ ] 10.2 Add performance monitoring
    - Add timing measurements for major operations (thermodynamic calculations, mass gap solving, integration)
    - Log performance metrics: calculation time, convergence iterations, memory allocation
    - Add @info logs with elapsed_ms for completed operations
    - _Requirements: 15.7_
  
  - [ ] 10.3 Implement log configuration support
    - Document how to set global log level using Logging.global_logger
    - Document how to enable DEBUG logging for specific modules using ENV["JULIA_DEBUG"]
    - Document how to log to files using LoggingExtras
    - Add examples to documentation
    - _Requirements: 15.8_
  
  - [ ]* 10.4 Write unit tests for logging
    - Test that appropriate log messages are generated at each level
    - Test that structured logging format is used
    - Test that large arrays are not logged (use summary statistics)
    - _Requirements: 14.7_

- [ ] 11. Create documentation
  - [ ] 11.1 Write API documentation
    - Document all Abstract_Interface methods with DocStrings
    - Document all base utility functions with parameter descriptions and return types
    - Document exception types and when they are thrown
    - Create `docs/api/data_contracts.md` documenting parameter formats and conversion rules
    - _Requirements: 10.1_
  
  - [ ] 11.2 Write developer guides
    - Create `docs/guides/error_handling.md` explaining exception hierarchy and fallback strategies
    - Create `docs/guides/configuration.md` explaining log configuration and parameter management
    - Create `docs/guides/migration_guide.md` for transitioning from legacy API to new architecture
    - _Requirements: 10.4_
  
  - [ ] 11.3 Write examples
    - Create example showing how to implement a new concrete model
    - Create example showing how to use factory pattern to create models
    - Create example showing how to switch between PNJL and rPNJL models
    - Create example showing parameter format conversions
    - _Requirements: 10.2, 10.3_
  
  - [ ] 11.4 Document correctness properties
    - Document all 15 correctness properties in design document
    - Explain what each property validates
    - Link properties to their corresponding property-based tests
    - _Requirements: 10.5_

- [ ] 12. Performance optimization and validation
  - [ ] 12.1 Add performance hints
    - Add @inline hints to performance-critical functions (dispersion_relation, distribution functions, energy calculations)
    - Ensure type stability throughout critical paths
    - Use @inbounds where safe for array access
    - _Requirements: 7.2_
  
  - [ ] 12.2 Implement closure-based optimization
    - Return closures for repeated calculations to enable JIT optimization
    - Cache expensive computations where appropriate
    - _Requirements: 7.3_
  
  - [ ] 12.3 Create performance benchmark suite
    - Create `tests/perf/models/benchmark_pnjl.jl` comparing refactored vs original PNJL
    - Create `tests/perf/models/benchmark_rpnjl.jl` for rPNJL performance
    - Measure critical path performance: thermodynamic calculations, mass gap solving, integration
    - Generate performance report comparing to baseline
    - _Requirements: 7.1, 7.5_
  
  - [ ]* 12.4 Run performance benchmarks
    - Execute benchmark suite and collect timing data
    - Compare refactored PNJL performance to original implementation
    - Verify performance degradation < 5%
    - Generate performance report
    - _Requirements: 7.5_

- [ ] 13. Integration testing and final validation
  - [ ]* 13.1 Write integration tests
    - Create `tests/integration/test_full_workflow.jl` testing complete thermodynamic calculation workflows
    - Test end-to-end: model creation → parameter setting → calculation → result validation
    - Test model switching: create PNJL model, run calculation, create rPNJL model, run calculation
    - Test with various parameter combinations and edge cases
    - _Requirements: 8.4_
  
  - [ ]* 13.2 Write edge case tests
    - Create `tests/edge_cases/test_zero_temperature.jl` for T=0 edge cases
    - Create `tests/edge_cases/test_extreme_mu.jl` for extreme chemical potential values
    - Create `tests/edge_cases/test_numerical_limits.jl` for numerical boundary conditions
    - _Requirements: 8.6_
  
  - [ ]* 13.3 Validate test coverage
    - Run coverage analysis on test suite
    - Verify core functionality has ≥80% coverage
    - Verify base utilities have 100% coverage
    - Generate coverage report
    - _Requirements: 8.1, 8.2_

- [ ] 14. Final checkpoint - Complete validation
  - Ensure all unit tests pass
  - Ensure all 15 property-based tests pass (Properties 1-15)
  - Ensure all regression tests pass with <1e-12 error
  - Ensure all integration tests pass
  - Ensure performance benchmarks show <5% degradation
  - Ensure test coverage ≥80% for core functionality
  - Ensure all documentation is complete
  - Ask the user if questions arise

## Notes

- Tasks marked with `*` are optional and can be skipped for faster MVP
- Each task references specific requirements for traceability
- Property-based tests validate all 15 correctness properties from the design document
- Regression tests ensure numerical equivalence with original implementation
- Checkpoints ensure incremental validation at key milestones
- Implementation follows phased approach: Foundation → PNJL → rPNJL → Integration
