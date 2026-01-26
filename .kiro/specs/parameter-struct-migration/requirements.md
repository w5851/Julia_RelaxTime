# Requirements Document: Parameter Struct Migration

## Introduction

This specification defines the requirements for migrating the PNJL (Polyakov-Nambu-Jona-Lasinio) model Julia codebase from mixed NamedTuple/struct parameter passing to a unified struct-based approach. The migration focuses on Phase A: full struct-ification of the RelaxationTime module chain while maintaining backward compatibility with existing NamedTuple-based code.

The current codebase suffers from type inconsistency issues when the same struct is defined in multiple modules, difficulty maintaining parameter validation and defaults, and verbose function signatures. This migration will establish `QuarkParams` and `ThermoParams` as the standard parameter representation across the project while preserving all existing functionality.

## Glossary

- **RelaxationTime_Module**: The main module (`src/relaxtime/RelaxationTime.jl`) that calculates quark relaxation times based on average scattering rates
- **QuarkParams**: Project-level struct encapsulating quark masses `m=(u,d,s)` and chemical potentials `μ=(u,d,s)`
- **ThermoParams**: Project-level struct encapsulating thermodynamic parameters `T, Φ, Φbar, ξ`
- **Dual_Interface**: API design pattern where public functions accept both `Union{NamedTuple,QuarkParams}` and `Union{NamedTuple,ThermoParams}`
- **Internal_Normalization**: Pattern using `_nt_quark(q)` and `_nt_thermo(t)` helper functions at function boundaries to convert inputs to consistent internal representation
- **AverageScatteringRate_Module**: Module (`src/relaxtime/AverageScatteringRate.jl`) that computes average scattering rates for quark-quark interactions
- **TotalCrossSection_Module**: Module (`src/relaxtime/TotalCrossSection.jl`) that calculates total scattering cross-sections by integrating differential cross-sections
- **ScatteringAmplitude_Module**: Module (`src/relaxtime/ScatteringAmplitude.jl`) that computes scattering amplitudes for quark interactions
- **DifferentialCrossSection_Module**: Module (`src/relaxtime/DifferentialCrossSection.jl`) that calculates differential scattering cross-sections
- **CrossSectionCache**: Struct that stores precomputed cross-section values with PCHIP interpolation for performance optimization
- **PNJL_Cutoff**: Momentum cutoff parameter Λ (approximately 3 fm⁻¹) used in PNJL model calculations
- **TotalPropagator_Module**: Module (`src/relaxtime/TotalPropagator.jl`) that calculates total propagators for meson exchange in scattering processes
- **ParticleSymbols_Module**: Utility module (`src/utils/ParticleSymbols.jl`) that provides particle symbol parsing and parameter lookup functions
- **TransportCoefficients_Module**: Module (`src/relaxtime/TransportCoefficients.jl`) that calculates transport coefficients (already has struct support via TransportRequest)

## Requirements

### Requirement 1: Dual Interface Support for RelaxationTime Module

**User Story:** As a developer using the RelaxationTime module, I want to pass parameters as either structs or NamedTuples, so that I can migrate gradually without breaking existing code.

#### Acceptance Criteria

1. WHEN a user calls `relaxation_times` with `QuarkParams` and `ThermoParams` structs, THE RelaxationTime_Module SHALL compute relaxation times correctly
2. WHEN a user calls `relaxation_times` with NamedTuple parameters, THE RelaxationTime_Module SHALL compute identical results to struct inputs
3. WHEN a user calls `compute_average_rates` with struct parameters, THE RelaxationTime_Module SHALL accept and process them correctly
4. WHEN a user calls `relaxation_rates` with struct-derived densities, THE RelaxationTime_Module SHALL compute rates correctly
5. THE RelaxationTime_Module SHALL use `_nt_quark` and `_nt_thermo` helper functions to normalize inputs at function boundaries

### Requirement 2: Dual Interface Support for AverageScatteringRate Module

**User Story:** As a developer calling average scattering rate functions, I want to pass parameters as either structs or NamedTuples, so that my code remains compatible during the migration.

#### Acceptance Criteria

1. WHEN a user calls `average_scattering_rate` with `QuarkParams` and `ThermoParams` structs, THE AverageScatteringRate_Module SHALL compute scattering rates correctly
2. WHEN a user calls `average_scattering_rate` with NamedTuple parameters, THE AverageScatteringRate_Module SHALL produce identical results to struct inputs
3. THE AverageScatteringRate_Module SHALL normalize struct inputs to NamedTuples at function entry using helper functions
4. WHEN `get_mass` or `get_mu` helper functions receive struct parameters, THE AverageScatteringRate_Module SHALL extract values correctly
5. WHEN `build_w0cdf_pchip_cache` receives struct parameters, THE AverageScatteringRate_Module SHALL construct CrossSectionCache correctly

### Requirement 3: Dual Interface Support for TotalCrossSection Module

**User Story:** As a developer computing cross-sections, I want to pass parameters as either structs or NamedTuples, so that I can use the new struct interface without modifying existing code.

#### Acceptance Criteria

1. WHEN a user calls `total_cross_section` with `QuarkParams` and `ThermoParams` structs, THE TotalCrossSection_Module SHALL compute cross-sections correctly
2. WHEN a user calls `total_cross_section` with NamedTuple parameters, THE TotalCrossSection_Module SHALL produce identical results to struct inputs
3. THE TotalCrossSection_Module SHALL normalize struct inputs at function boundaries
4. WHEN `calculate_all_total_cross_sections` receives struct parameters, THE TotalCrossSection_Module SHALL process all scattering processes correctly
5. WHEN `scan_s_dependence` receives struct parameters, THE TotalCrossSection_Module SHALL generate correct s-dependent cross-section data

### Requirement 4: Dual Interface Support for ScatteringAmplitude Module

**User Story:** As a developer working with scattering amplitudes, I want to pass parameters as either structs or NamedTuples, so that amplitude calculations work with both parameter formats.

#### Acceptance Criteria

1. WHEN a user calls `scattering_amplitude_squared` with `QuarkParams` and `ThermoParams` structs, THE ScatteringAmplitude_Module SHALL compute amplitudes correctly
2. WHEN a user calls `scattering_amplitude_squared` with NamedTuple parameters, THE ScatteringAmplitude_Module SHALL produce identical results to struct inputs
3. THE ScatteringAmplitude_Module SHALL normalize struct inputs at function entry
4. WHEN amplitude functions access quark masses or chemical potentials from structs, THE ScatteringAmplitude_Module SHALL extract values correctly
5. WHEN amplitude functions receive K_coeffs alongside struct parameters, THE ScatteringAmplitude_Module SHALL process coupling coefficients correctly

### Requirement 5: Dual Interface Support for DifferentialCrossSection Module

**User Story:** As a developer computing differential cross-sections, I want to pass parameters as either structs or NamedTuples, so that differential cross-section calculations support both formats.

#### Acceptance Criteria

1. WHEN a user calls `differential_cross_section` with `QuarkParams` and `ThermoParams` structs, THE DifferentialCrossSection_Module SHALL compute differential cross-sections correctly
2. WHEN a user calls `differential_cross_section` with NamedTuple parameters, THE DifferentialCrossSection_Module SHALL produce identical results to struct inputs
3. THE DifferentialCrossSection_Module SHALL normalize struct inputs at function boundaries
4. WHEN differential cross-section functions call scattering amplitude functions with structs, THE DifferentialCrossSection_Module SHALL pass parameters correctly
5. WHEN kinematic calculations use struct parameters, THE DifferentialCrossSection_Module SHALL compute Mandelstam variables correctly

### Requirement 6: Dual Interface Support for TotalPropagator Module

**User Story:** As a developer computing propagators for meson exchange, I want to pass parameters as either structs or NamedTuples, so that propagator calculations support both formats.

#### Acceptance Criteria

1. WHEN a user calls `total_propagator_simple` with `QuarkParams` and `ThermoParams` structs, THE TotalPropagator_Module SHALL compute propagators correctly
2. WHEN a user calls `total_propagator_mixed` with struct parameters, THE TotalPropagator_Module SHALL produce identical results to NamedTuple inputs
3. WHEN a user calls `calculate_cms_momentum` with struct parameters, THE TotalPropagator_Module SHALL compute quality center-of-mass momentum correctly
4. WHEN a user calls `calculate_all_propagators` with struct parameters, THE TotalPropagator_Module SHALL compute all relevant propagators correctly
5. THE TotalPropagator_Module SHALL normalize struct inputs at function boundaries

### Requirement 7: Dual Interface Support for ParticleSymbols Module

**User Story:** As a developer using particle symbol utilities, I want parameter lookup functions to accept both structs and NamedTuples, so that symbol utilities work seamlessly with the new struct interface.

#### Acceptance Criteria

1. WHEN a user calls `get_mass` with `QuarkParams` struct, THE ParticleSymbols_Module SHALL extract particle masses correctly
2. WHEN a user calls `get_chemical_potential` with `QuarkParams` struct, THE ParticleSymbols_Module SHALL extract chemical potentials correctly
3. WHEN a user calls `get_quark_masses_for_process` with struct parameters, THE ParticleSymbols_Module SHALL return correct mass tuples for all four particles
4. THE ParticleSymbols_Module SHALL normalize struct inputs using helper functions
5. WHEN struct and NamedTuple inputs are used interchangeably, THE ParticleSymbols_Module SHALL produce identical results

### Requirement 8: Backward Compatibility Preservation

**User Story:** As a developer with existing scripts and tests using NamedTuples, I want all my code to continue working without modification, so that the migration does not break my workflow.

#### Acceptance Criteria

1. WHEN existing test suites run without modification, THE System SHALL pass all tests that previously passed
2. WHEN existing scripts call functions with NamedTuple parameters, THE System SHALL produce identical results to pre-migration behavior
3. WHEN users mix struct and NamedTuple parameters in the same workflow, THE System SHALL handle conversions transparently
4. THE System SHALL NOT introduce any breaking API changes to public function signatures
5. WHEN documentation examples use NamedTuples, THE System SHALL continue to support those usage patterns

### Requirement 9: Internal Consistency and Validation

**User Story:** As a developer, I want struct parameters to be validated on construction and used consistently internally, so that I catch errors early and maintain code quality.

#### Acceptance Criteria

1. WHEN a QuarkParams struct is constructed with missing required fields, THE System SHALL raise a descriptive error
2. WHEN a ThermoParams struct is constructed with missing required fields, THE System SHALL raise a descriptive error
3. WHEN internal functions receive normalized parameters, THE System SHALL use a consistent representation throughout the call chain
4. WHEN `ensure_quark_params_has_A` processes struct inputs, THE System SHALL compute and attach the A field correctly
5. THE System SHALL provide `as_namedtuple` conversion functions for backward compatibility with legacy interfaces

### Requirement 10: Performance Preservation

**User Story:** As a developer running physics simulations, I want the struct migration to maintain or improve performance, so that my computational workflows remain efficient.

#### Acceptance Criteria

1. WHEN relaxation time calculations use struct parameters, THE System SHALL execute within 5% of NamedTuple performance
2. WHEN average scattering rate calculations use struct parameters, THE System SHALL execute within 5% of NamedTuple performance
3. WHEN cross-section calculations use struct parameters, THE System SHALL execute within 5% of NamedTuple performance
4. WHEN helper functions normalize inputs, THE System SHALL use inlined functions to minimize overhead
5. WHEN CrossSectionCache operations use struct parameters, THE System SHALL maintain cache hit rates and interpolation performance

### Requirement 11: Test Coverage for Struct Equivalence

**User Story:** As a developer ensuring correctness, I want comprehensive tests verifying that struct inputs produce identical results to NamedTuple inputs, so that I can trust the migration.

#### Acceptance Criteria

1. WHEN unit tests compare struct vs NamedTuple inputs for `relaxation_times`, THE System SHALL verify numerical equivalence within tolerance
2. WHEN unit tests compare struct vs NamedTuple inputs for `average_scattering_rate`, THE System SHALL verify numerical equivalence within tolerance
3. WHEN unit tests compare struct vs NamedTuple inputs for `total_cross_section`, THE System SHALL verify numerical equivalence within tolerance
4. WHEN unit tests compare struct vs NamedTuple inputs for `scattering_amplitude_squared`, THE System SHALL verify numerical equivalence within tolerance
5. WHEN unit tests compare struct vs NamedTuple inputs for `differential_cross_section`, THE System SHALL verify numerical equivalence within tolerance
6. WHEN unit tests compare struct vs NamedTuple inputs for `total_propagator_simple`, THE System SHALL verify numerical equivalence within tolerance
7. WHEN unit tests compare struct vs NamedTuple inputs for `get_mass` and `get_chemical_potential`, THE System SHALL verify correct value extraction

### Requirement 12: Documentation and Usage Guidance

**User Story:** As a developer learning the codebase, I want clear documentation showing how to use struct parameters, so that I can adopt best practices in my code.

#### Acceptance Criteria

1. WHEN a developer reads function docstrings, THE System SHALL document both struct and NamedTuple parameter formats
2. WHEN a developer reads module documentation, THE System SHALL explain the dual interface pattern and internal normalization strategy
3. WHEN a developer reads migration guides, THE System SHALL provide examples of converting NamedTuple code to struct code
4. WHEN a developer reads API documentation, THE System SHALL show struct usage as the primary recommended pattern
5. WHEN a developer encounters helper functions like `_nt_quark`, THE System SHALL document their purpose and usage
