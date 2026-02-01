# Design Document: QCD Model Architecture Refactoring

## Overview

This design document specifies the technical architecture for refactoring the QCD model library from a monolithic PNJL implementation to a polymorphic, extensible system supporting multiple QCD models (PNJL, rPNJL, magnetic PNJL, etc.).

### Design Goals

1. **Polymorphic Architecture**: Leverage Julia's multiple dispatch to create a flexible type hierarchy that allows easy addition of new models
2. **Code Reuse**: Extract common numerical utilities into shared base modules to eliminate duplication
3. **Backward Compatibility**: Maintain existing API surface so current scripts continue working without modification
4. **Performance**: Ensure refactored code maintains ≤5% performance degradation compared to original implementation
5. **Extensibility**: Enable researchers to add new models by implementing a small set of required methods

### Key Design Principles

- **Interface Segregation**: Define minimal abstract interfaces that models must implement
- **Dependency Inversion**: High-level thermodynamic calculations depend on abstract interfaces, not concrete implementations
- **Open/Closed**: System is open for extension (new models) but closed for modification (core infrastructure)
- **Type Stability**: Maintain type stability throughout for Julia compiler optimization
- **Fail-Safe Defaults**: Provide safe numerical operations with graceful degradation

### Architecture Style

The design follows a **layered architecture** with clear separation of concerns:


```
┌─────────────────────────────────────────────────────────┐
│           User Scripts / Legacy API                     │
│         (Compatibility Layer - Optional)                │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│              High-Level API / Factory                   │
│         (create_model, parameter conversion)            │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│            Abstract Model Interface                     │
│    (AbstractQCDModel, domain-specific abstracts)        │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│           Concrete Model Implementations                │
│         (PNJLModel, RPNJLModel, etc.)                   │
└─────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────┐
│              Base Utilities Layer                       │
│  (Integrals, Polyakov, SafeMath, Distributions)         │
└─────────────────────────────────────────────────────────┘
```


## Architecture

### Type Hierarchy

The system uses Julia's abstract type system to define a hierarchy of model types:

```julia
abstract type AbstractQCDModel end

# Domain-specific abstract types
abstract type AbstractIsotropicModel <: AbstractQCDModel end
abstract type AbstractAnisotropicModel <: AbstractQCDModel end
abstract type AbstractMagneticModel <: AbstractQCDModel end

# Concrete implementations
struct PNJLModel <: AbstractIsotropicModel
    params::QuarkParams
    thermo_params::ThermoParams
    # ... additional fields
end

struct RPNJLModel <: AbstractAnisotropicModel
    params::QuarkParams
    thermo_params::ThermoParams
    g1::Float64  # Eight-quark coupling
    g2::Float64  # Eight-quark coupling
    kappa::Float64  # Vandermonde parameter
    T0::Float64  # Reference temperature
    # ... additional fields
end
```

### Module Organization

The codebase is organized into a clear directory structure:

```
src/
├── core/
│   ├── interface.jl          # Abstract types and required method signatures
│   └── exceptions.jl          # Custom exception hierarchy
├── models/
│   ├── base/                  # Shared utilities
│   │   ├── integrals.jl       # Numerical integration
│   │   ├── polyakov.jl        # Polyakov loop calculations
│   │   ├── safe_math.jl       # Safe numerical operations
│   │   ├── distributions.jl   # Quark/antiquark distributions
│   │   ├── validation.jl      # Parameter validation
│   │   └── mass_chiral.jl     # Mass and chiral condensate
│   ├── isotropic/
│   │   ├── types.jl           # AbstractIsotropicModel
│   │   └── pnjl.jl            # PNJLModel implementation
│   ├── anisotropic/
│   │   ├── types.jl           # AbstractAnisotropicModel
│   │   └── rpnjl.jl           # RPNJLModel implementation
│   ├── magnetic/
│   │   └── types.jl           # AbstractMagneticModel (future)
│   └── factory.jl             # Model creation factory
└── compatibility/
    └── legacy_api.jl          # Backward compatibility layer
```


### Multi-Dispatch Strategy

Julia's multiple dispatch enables elegant polymorphism without traditional OOP patterns:

```julia
# Generic interface - works for all AbstractQCDModel subtypes
function calculate_thermodynamics(model::AbstractQCDModel, T::Float64, μ::Float64)
    # Dispatch to model-specific implementations
    vac = vacuum_contribution(model)
    therm = thermal_contribution(model, T, μ)
    poly = polyakov_potential(model, T)
    return vac + therm + poly
end

# Model-specific implementations
function thermal_contribution(model::PNJLModel, T::Float64, μ::Float64)
    # PNJL-specific calculation
end

function thermal_contribution(model::RPNJLModel, T::Float64, μ::Float64)
    # rPNJL-specific calculation with eight-quark terms
end
```

This approach provides:
- **Type-based dispatch**: Compiler selects correct method at compile time
- **Zero runtime overhead**: No virtual function calls or dynamic dispatch
- **Extensibility**: New models automatically work with existing high-level functions
- **Specialization**: Compiler generates optimized code for each concrete type


## Components and Interfaces

### Core Interface Layer

**File**: `src/core/interface.jl`

Defines the abstract type hierarchy and required method signatures:

```julia
"""
Abstract base type for all QCD models.
All concrete models must inherit from this type and implement required methods.
"""
abstract type AbstractQCDModel end

# Required methods that all models must implement
"""
Calculate vacuum contribution to thermodynamic potential.
Returns: Float64 - vacuum energy density
"""
function vacuum_contribution(model::AbstractQCDModel)
    error("vacuum_contribution not implemented for $(typeof(model))")
end

"""
Calculate thermal contribution to thermodynamic potential.
Args:
    model: QCD model instance
    T: Temperature (MeV)
    μ: Chemical potential (MeV)
Returns: Float64 - thermal energy density
"""
function thermal_contribution(model::AbstractQCDModel, T::Float64, μ::Float64)
    error("thermal_contribution not implemented for $(typeof(model))")
end

"""
Calculate Polyakov loop potential.
Args:
    model: QCD model instance
    T: Temperature (MeV)
Returns: Float64 - Polyakov potential
"""
function polyakov_potential(model::AbstractQCDModel, T::Float64)
    error("polyakov_potential not implemented for $(typeof(model))")
end

"""
Calculate quark dispersion relation.
Args:
    model: QCD model instance
    p: Momentum (MeV)
    M: Constituent quark mass (MeV)
Returns: Float64 - quark energy
"""
function dispersion_relation(model::AbstractQCDModel, p::Float64, M::Float64)
    error("dispersion_relation not implemented for $(typeof(model))")
end

"""
Calculate constituent quark mass vector for all flavors.
Args:
    model: QCD model instance
    phi: Polyakov loop expectation value
    T: Temperature (MeV)
    μ: Chemical potential (MeV)
Returns: Vector{Float64} - mass for each flavor
"""
function calculate_mass_vec(model::AbstractQCDModel, phi::Float64, T::Float64, μ::Float64)
    error("calculate_mass_vec not implemented for $(typeof(model))")
end

"""
Calculate chiral condensate for all flavors.
Args:
    model: QCD model instance
    M_vec: Constituent quark masses (MeV)
    T: Temperature (MeV)
    μ: Chemical potential (MeV)
Returns: Vector{Float64} - chiral condensate for each flavor
"""
function calculate_chiral(model::AbstractQCDModel, M_vec::Vector{Float64}, T::Float64, μ::Float64)
    error("calculate_chiral not implemented for $(typeof(model))")
end
```


### Exception Hierarchy

**File**: `src/core/exceptions.jl`

Custom exception types for clear error handling:

```julia
"""Base exception for all QCD model errors"""
abstract type QCDModelException <: Exception end

"""Exception for invalid model parameters"""
struct ModelParameterError <: QCDModelException
    parameter::String
    value::Any
    constraint::String
    message::String
end

function Base.showerror(io::IO, e::ModelParameterError)
    print(io, "ModelParameterError: Invalid parameter '$(e.parameter)' = $(e.value)\n")
    print(io, "Constraint: $(e.constraint)\n")
    print(io, "Suggestion: $(e.message)")
end

"""Exception for numerical convergence failures"""
struct NumericalConvergenceError <: QCDModelException
    operation::String
    iterations::Int
    residual::Float64
    message::String
end

"""Exception for physical constraint violations"""
struct PhysicalConstraintError <: QCDModelException
    constraint::String
    value::Any
    message::String
end

"""Exception for model configuration errors"""
struct ModelConfigurationError <: QCDModelException
    model_type::String
    message::String
end

"""Exception for numerical instability"""
struct NumericalInstabilityError <: QCDModelException
    operation::String
    input::Any
    message::String
end
```


### Base Utilities Layer

#### Numerical Integration Module

**File**: `src/models/base/integrals.jl`

Provides reusable numerical integration with callback-based energy functions:

```julia
"""
Generic momentum integration for thermodynamic quantities.
Uses Gauss-Legendre quadrature with adaptive node count.

Args:
    energy_func: Function p -> E(p) that computes energy at momentum p
    T: Temperature (MeV)
    μ: Chemical potential (MeV)
    n_nodes: Number of quadrature nodes (default: 64)
Returns:
    Float64 - integrated value
"""
function integrate_momentum(energy_func::Function, T::Float64, μ::Float64; n_nodes::Int=64)
    # Implementation with fallback strategy
    try
        result = gauss_legendre_integrate(energy_func, 0.0, cutoff, n_nodes)
        return result
    catch e
        @warn "Integration failed with $n_nodes nodes, trying $(2*n_nodes)" exception=e
        return gauss_legendre_integrate(energy_func, 0.0, cutoff, 2*n_nodes)
    end
end

"""
Vacuum integral with UV cutoff regularization.
"""
function vacuum_integral(energy_func::Function, cutoff::Float64; n_nodes::Int=64)
    # Gauss-Legendre integration from 0 to cutoff
end
```


#### Polyakov Loop Module

**File**: `src/models/base/polyakov.jl`

Shared Polyakov loop calculations:

```julia
"""
Calculate Polyakov loop logarithm: ln(Φ) + ln(Φ*)/3
Used in thermal distribution functions.

Args:
    phi: Polyakov loop expectation value [0, 1]
    phi_bar: Conjugate Polyakov loop expectation value [0, 1]
Returns:
    Float64 - logarithmic combination
"""
function polyakov_log_combination(phi::Float64, phi_bar::Float64)
    # Validate inputs
    if !(0 <= phi <= 1) || !(0 <= phi_bar <= 1)
        throw(PhysicalConstraintError(
            "Polyakov loop range",
            (phi, phi_bar),
            "Polyakov loops must be in [0, 1]"
        ))
    end
    
    # Use safe_log to handle edge cases
    return safe_log(phi) + safe_log(phi_bar) / 3
end

"""
Standard Polyakov potential (polynomial form).
Used in PNJL model.
"""
function standard_polyakov_potential(phi::Float64, T::Float64, params)
    # Polynomial potential implementation
end

"""
Vandermonde-modified Polyakov potential.
Used in rPNJL model with temperature-dependent b2(T) term.

Args:
    phi: Polyakov loop expectation value
    T: Temperature (MeV)
    kappa: Vandermonde parameter
    T0: Reference temperature (MeV)
Returns:
    Float64 - Polyakov potential with Vandermonde terms
"""
function vandermonde_polyakov_potential(phi::Float64, T::Float64, kappa::Float64, T0::Float64)
    # Calculate b2(T) temperature dependence
    b2_T = calculate_b2_temperature(T, T0, kappa)
    
    # Compute Vandermonde terms (equations 3.27-3.29)
    vander_term = safe_vandermonde(phi, b2_T)
    
    # Combine with standard potential
    return standard_polyakov_potential(phi, T, params) + vander_term
end
```


#### Safe Math Module

**File**: `src/models/base/safe_math.jl`

Numerical operations with graceful error handling:

```julia
"""
Safe logarithm that handles edge cases.

Args:
    x: Input value
    fallback: Value to return for invalid inputs (default: -Inf)
Returns:
    Float64 - log(x) or fallback value
"""
function safe_log(x::Float64; fallback::Float64=-Inf)
    if x <= 0.0
        if x == 0.0
            return fallback
        else
            throw(NumericalInstabilityError(
                "logarithm",
                x,
                "Cannot take logarithm of negative number: $x"
            ))
        end
    end
    return log(x)
end

"""
Safe Vandermonde term calculation with negative discriminant handling.

Args:
    phi: Polyakov loop value
    b2: Temperature-dependent parameter
Returns:
    Float64 - Vandermonde term value
    
Note: Logs warning if discriminant is negative but continues calculation.
"""
function safe_vandermonde(phi::Float64, b2::Float64)
    discriminant = compute_discriminant(phi, b2)
    
    if discriminant < 0
        @warn "Negative Vandermonde discriminant" phi=phi b2=b2 discriminant=discriminant
        # Continue with calculation - physical interpretation may still be valid
    end
    
    return calculate_vandermonde_term(phi, b2, discriminant)
end

"""
Safe division with zero check.
"""
function safe_divide(numerator::Float64, denominator::Float64; fallback::Float64=0.0)
    if abs(denominator) < eps(Float64)
        @warn "Division by near-zero denominator" numerator=numerator denominator=denominator
        return fallback
    end
    return numerator / denominator
end
```


#### Distribution Functions Module

**File**: `src/models/base/distributions.jl`

Fermi-Dirac distributions for quarks and antiquarks:

```julia
"""
Quark Fermi-Dirac distribution with Polyakov loop modification.

Args:
    E: Quark energy (MeV)
    T: Temperature (MeV)
    μ: Chemical potential (MeV)
    phi: Polyakov loop expectation value
Returns:
    Float64 - distribution function value
"""
function quark_distribution(E::Float64, T::Float64, μ::Float64, phi::Float64)
    if T <= 0.0
        # Zero temperature limit
        return E < μ ? 1.0 : 0.0
    end
    
    # Polyakov-modified distribution
    exponent = (E - μ) / T
    return phi / (1.0 + phi * exp(exponent))
end

"""
Antiquark Fermi-Dirac distribution with Polyakov loop modification.

Args:
    E: Antiquark energy (MeV)
    T: Temperature (MeV)
    μ: Chemical potential (MeV)
    phi_bar: Conjugate Polyakov loop expectation value
Returns:
    Float64 - distribution function value
"""
function antiquark_distribution(E::Float64, T::Float64, μ::Float64, phi_bar::Float64)
    if T <= 0.0
        return 0.0  # No antiquarks at T=0
    end
    
    exponent = (E + μ) / T
    return phi_bar / (1.0 + phi_bar * exp(exponent))
end
```


#### Parameter Validation Module

**File**: `src/models/base/validation.jl`

Input validation with clear constraint checking:

```julia
"""
Validate model parameters against physical constraints.

Args:
    params: Parameter structure or dictionary
    constraints: Dictionary of parameter_name => (min, max, description)
Throws:
    ModelParameterError if any constraint is violated
"""
function validate_parameters(params, constraints::Dict)
    for (param_name, (min_val, max_val, description)) in constraints
        value = get_parameter_value(params, param_name)
        
        if value < min_val || value > max_val
            throw(ModelParameterError(
                string(param_name),
                value,
                "$min_val <= $param_name <= $max_val",
                "Parameter $param_name must satisfy: $description"
            ))
        end
    end
end

"""
Standard constraints for QCD model parameters.
"""
const STANDARD_CONSTRAINTS = Dict(
    :T => (0.0, 500.0, "Temperature in MeV, physical range 0-500 MeV"),
    :mu => (-500.0, 500.0, "Chemical potential in MeV"),
    :phi => (0.0, 1.0, "Polyakov loop expectation value"),
    :M => (0.0, 1000.0, "Constituent quark mass in MeV"),
)

"""
Validate physical constraints during calculation.
"""
function validate_physical_constraint(name::String, value::Float64, condition::Bool, message::String)
    if !condition
        throw(PhysicalConstraintError(name, value, message))
    end
end
```


#### Mass and Chiral Condensate Module

**File**: `src/models/base/mass_chiral.jl`

Shared calculations for mass gap equations and chiral condensates:

```julia
"""
Solve mass gap equation using iterative method.

Args:
    model: QCD model instance
    initial_mass: Initial guess for mass (MeV)
    T: Temperature (MeV)
    μ: Chemical potential (MeV)
    phi: Polyakov loop value
    max_iter: Maximum iterations (default: 100)
    tol: Convergence tolerance (default: 1e-6)
Returns:
    Float64 - converged constituent mass
Throws:
    NumericalConvergenceError if solver fails to converge
"""
function solve_mass_gap(model::AbstractQCDModel, initial_mass::Float64, T::Float64, μ::Float64, phi::Float64;
                        max_iter::Int=100, tol::Float64=1e-6)
    M = initial_mass
    
    for iter in 1:max_iter
        M_new = calculate_mass_iteration(model, M, T, μ, phi)
        residual = abs(M_new - M)
        
        if residual < tol
            @debug "Mass gap converged" iterations=iter residual=residual
            return M_new
        end
        
        M = M_new
    end
    
    # Convergence failed - try fallback strategies
    @warn "Mass gap equation did not converge, attempting fallback" max_iter=max_iter
    
    # Strategy 1: Relax tolerance
    return solve_mass_gap(model, initial_mass, T, μ, phi; max_iter=max_iter*2, tol=tol*10)
end

"""
Calculate chiral condensate from constituent mass.

Args:
    model: QCD model instance
    M: Constituent quark mass (MeV)
    T: Temperature (MeV)
    μ: Chemical potential (MeV)
Returns:
    Float64 - chiral condensate value
"""
function calculate_chiral_condensate(model::AbstractQCDModel, M::Float64, T::Float64, μ::Float64)
    # Integration over momentum with distribution functions
    integrand = p -> chiral_integrand(model, p, M, T, μ)
    return integrate_momentum(integrand, T, μ)
end
```


### Concrete Model Implementations

#### PNJL Model

**File**: `src/models/isotropic/pnjl.jl`

Refactored PNJL model using new architecture:

```julia
"""
Polyakov-loop extended Nambu-Jona-Lasinio model.
Isotropic model with standard Polyakov potential.
"""
struct PNJLModel <: AbstractIsotropicModel
    params::QuarkParams
    thermo_params::ThermoParams
    cutoff::Float64
    
    # Constructor with validation
    function PNJLModel(params::QuarkParams, thermo_params::ThermoParams, cutoff::Float64)
        validate_parameters(params, STANDARD_CONSTRAINTS)
        validate_parameters(thermo_params, STANDARD_CONSTRAINTS)
        new(params, thermo_params, cutoff)
    end
end

# Implement required interface methods
function vacuum_contribution(model::PNJLModel)
    # Use base utilities for integration
    integrand = p -> vacuum_energy_density(p, model.params)
    return vacuum_integral(integrand, model.cutoff)
end

function thermal_contribution(model::PNJLModel, T::Float64, μ::Float64)
    # Solve for constituent mass
    phi = calculate_polyakov_loop(model, T)
    M_vec = calculate_mass_vec(model, phi, T, μ)
    
    # Integrate thermal contribution
    integrand = p -> thermal_energy_density(p, M_vec, T, μ, phi)
    return integrate_momentum(integrand, T, μ)
end

function polyakov_potential(model::PNJLModel, T::Float64)
    phi = calculate_polyakov_loop(model, T)
    return standard_polyakov_potential(phi, T, model.thermo_params)
end

function dispersion_relation(model::PNJLModel, p::Float64, M::Float64)
    return sqrt(p^2 + M^2)
end

function calculate_mass_vec(model::PNJLModel, phi::Float64, T::Float64, μ::Float64)
    # Solve mass gap equation for each flavor
    masses = Float64[]
    for flavor in 1:model.params.n_flavors
        M_initial = model.params.current_masses[flavor]
        M = solve_mass_gap(model, M_initial, T, μ, phi)
        push!(masses, M)
    end
    return masses
end

function calculate_chiral(model::PNJLModel, M_vec::Vector{Float64}, T::Float64, μ::Float64)
    # Calculate chiral condensate for each flavor
    return [calculate_chiral_condensate(model, M, T, μ) for M in M_vec]
end
```


#### rPNJL Model

**File**: `src/models/anisotropic/rpnjl.jl`

New rPNJL model with eight-quark terms and Vandermonde potential:

```julia
"""
Repulsive PNJL model with eight-quark interaction terms.
Anisotropic model with modified Polyakov potential including Vandermonde terms.
"""
struct RPNJLModel <: AbstractAnisotropicModel
    params::QuarkParams
    thermo_params::ThermoParams
    cutoff::Float64
    
    # rPNJL-specific parameters
    g1::Float64      # Eight-quark coupling constant (dimensionless)
    g2::Float64      # Eight-quark coupling constant (dimensionless)
    kappa::Float64   # Vandermonde parameter (dimensionless)
    T0::Float64      # Reference temperature (MeV)
    
    function RPNJLModel(params::QuarkParams, thermo_params::ThermoParams, cutoff::Float64,
                        g1::Float64, g2::Float64, kappa::Float64, T0::Float64)
        # Validate standard parameters
        validate_parameters(params, STANDARD_CONSTRAINTS)
        validate_parameters(thermo_params, STANDARD_CONSTRAINTS)
        
        # Validate rPNJL-specific parameters
        rpnjl_constraints = Dict(
            :g1 => (0.0, 100.0, "Eight-quark coupling g1"),
            :g2 => (0.0, 100.0, "Eight-quark coupling g2"),
            :kappa => (0.0, 10.0, "Vandermonde parameter"),
            :T0 => (50.0, 300.0, "Reference temperature in MeV")
        )
        validate_parameters((g1=g1, g2=g2, kappa=kappa, T0=T0), rpnjl_constraints)
        
        new(params, thermo_params, cutoff, g1, g2, kappa, T0)
    end
end

# Implement required interface methods with rPNJL-specific physics

function vacuum_contribution(model::RPNJLModel)
    # Same as PNJL - vacuum term unchanged
    integrand = p -> vacuum_energy_density(p, model.params)
    return vacuum_integral(integrand, model.cutoff)
end

function thermal_contribution(model::RPNJLModel, T::Float64, μ::Float64)
    phi = calculate_polyakov_loop(model, T)
    M_vec = calculate_mass_vec(model, phi, T, μ)
    
    # Include eight-quark term contribution (equation 3.31)
    eight_quark_term = calculate_eight_quark_contribution(model, M_vec, T, μ)
    
    integrand = p -> thermal_energy_density(p, M_vec, T, μ, phi)
    thermal = integrate_momentum(integrand, T, μ)
    
    return thermal + eight_quark_term
end

function polyakov_potential(model::RPNJLModel, T::Float64)
    phi = calculate_polyakov_loop(model, T)
    # Use Vandermonde-modified potential (equations 3.27-3.29)
    return vandermonde_polyakov_potential(phi, T, model.kappa, model.T0)
end

function calculate_mass_vec(model::RPNJLModel, phi::Float64, T::Float64, μ::Float64)
    # Mass gap equation includes eight-quark terms (equation 3.31)
    masses = Float64[]
    for flavor in 1:model.params.n_flavors
        M_initial = model.params.current_masses[flavor]
        M = solve_mass_gap_rpnjl(model, M_initial, T, μ, phi)
        push!(masses, M)
    end
    return masses
end

"""
Calculate eight-quark interaction contribution (equation 3.31).
This is the key difference from standard PNJL.
"""
function calculate_eight_quark_contribution(model::RPNJLModel, M_vec::Vector{Float64}, T::Float64, μ::Float64)
    # Eight-quark term: -g1 * σ^4 - g2 * σ^2 * Δ^2
    # where σ is chiral condensate and Δ is diquark condensate
    
    chiral_vec = calculate_chiral(model, M_vec, T, μ)
    σ = sum(chiral_vec)
    
    # For now, assume Δ = 0 (no diquark condensate)
    # Full implementation would solve diquark gap equation
    
    return -model.g1 * σ^4 - model.g2 * σ^2 * 0.0^2
end

"""
Solve mass gap equation with eight-quark terms.
Modified from standard PNJL to include repulsive interactions.
"""
function solve_mass_gap_rpnjl(model::RPNJLModel, initial_mass::Float64, T::Float64, μ::Float64, phi::Float64;
                               max_iter::Int=100, tol::Float64=1e-6)
    M = initial_mass
    
    for iter in 1:max_iter
        # Standard PNJL term
        M_pnjl = calculate_mass_iteration(model, M, T, μ, phi)
        
        # Eight-quark correction
        eight_quark_correction = calculate_eight_quark_mass_correction(model, M, T, μ)
        
        M_new = M_pnjl + eight_quark_correction
        residual = abs(M_new - M)
        
        if residual < tol
            @debug "rPNJL mass gap converged" iterations=iter residual=residual
            return M_new
        end
        
        M = M_new
    end
    
    throw(NumericalConvergenceError(
        "rPNJL mass gap equation",
        max_iter,
        abs(M - initial_mass),
        "Try adjusting initial guess or increasing max_iter"
    ))
end
```


### Factory Pattern

**File**: `src/models/factory.jl`

Simplified model creation:

```julia
"""
Create a QCD model instance by name.

Args:
    model_name: Symbol identifying the model (:PNJL, :rPNJL, etc.)
    kwargs: Model-specific parameters
Returns:
    AbstractQCDModel - concrete model instance
Throws:
    ModelConfigurationError if model_name is unknown

Example:
    model = create_model(:PNJL, params=quark_params, thermo_params=thermo_params, cutoff=650.0)
    model = create_model(:rPNJL, params=quark_params, g1=10.0, g2=5.0, kappa=0.5, T0=270.0)
"""
function create_model(model_name::Symbol; kwargs...)
    if model_name == :PNJL
        return create_pnjl_model(; kwargs...)
    elseif model_name == :rPNJL
        return create_rpnjl_model(; kwargs...)
    elseif model_name == :MagneticPNJL
        return create_magnetic_pnjl_model(; kwargs...)
    else
        throw(ModelConfigurationError(
            string(model_name),
            "Unknown model type. Available models: :PNJL, :rPNJL, :MagneticPNJL"
        ))
    end
end

function create_pnjl_model(; params::QuarkParams, thermo_params::ThermoParams, cutoff::Float64=650.0)
    return PNJLModel(params, thermo_params, cutoff)
end

function create_rpnjl_model(; params::QuarkParams, thermo_params::ThermoParams, cutoff::Float64=650.0,
                             g1::Float64, g2::Float64, kappa::Float64, T0::Float64)
    return RPNJLModel(params, thermo_params, cutoff, g1, g2, kappa, T0)
end
```


### Backward Compatibility Layer

**File**: `src/compatibility/legacy_api.jl`

Maintains existing API for current scripts:

```julia
"""
Legacy API wrapper that translates old function calls to new architecture.
Ensures existing scripts continue working without modification.
"""

# Legacy function: calculate thermodynamics with old parameter format
function calculate_thermodynamics_legacy(T::Float64, μ::Float64, params_dict::Dict)
    # Convert legacy parameters to new format
    params = convert_legacy_params(params_dict)
    thermo_params = convert_legacy_thermo_params(params_dict)
    
    # Create model using factory
    model = create_model(:PNJL, params=params, thermo_params=thermo_params)
    
    # Call new API
    return calculate_thermodynamics(model, T, μ)
end

# Legacy function: direct PNJL calculation
function pnjl_thermodynamics(T::Float64, μ::Float64)
    # Use default parameters from legacy code
    params = get_default_pnjl_params()
    model = create_model(:PNJL, params=params.quark, thermo_params=params.thermo)
    return calculate_thermodynamics(model, T, μ)
end

"""
Convert legacy parameter dictionary to new QuarkParams struct.
Handles both old and new parameter names.
"""
function convert_legacy_params(params_dict::Dict)
    # Map old parameter names to new ones
    param_mapping = Dict(
        "m_current" => :current_masses,
        "G" => :coupling_G,
        "K" => :coupling_K,
        # ... additional mappings
    )
    
    converted = Dict{Symbol, Any}()
    for (old_name, new_name) in param_mapping
        if haskey(params_dict, old_name)
            converted[new_name] = params_dict[old_name]
        end
    end
    
    return QuarkParams(; converted...)
end

# Export legacy functions to maintain API compatibility
export calculate_thermodynamics_legacy, pnjl_thermodynamics
```


## Data Models

### Parameter Structures

```julia
"""
Quark sector parameters.
Immutable for thread safety.
"""
struct QuarkParams
    n_flavors::Int                    # Number of quark flavors (typically 2 or 3)
    current_masses::Vector{Float64}   # Current quark masses (MeV)
    coupling_G::Float64               # Four-quark coupling (GeV^-2)
    coupling_K::Float64               # Six-quark coupling (GeV^-5)
    cutoff::Float64                   # UV cutoff (MeV)
    
    # Constructor with validation
    function QuarkParams(n_flavors::Int, current_masses::Vector{Float64}, 
                        coupling_G::Float64, coupling_K::Float64, cutoff::Float64)
        @assert n_flavors > 0 "Number of flavors must be positive"
        @assert length(current_masses) == n_flavors "Mass vector length must match n_flavors"
        @assert coupling_G > 0 "Coupling G must be positive"
        @assert cutoff > 0 "Cutoff must be positive"
        new(n_flavors, current_masses, coupling_G, coupling_K, cutoff)
    end
end

"""
Thermodynamic parameters.
Immutable for thread safety.
"""
struct ThermoParams
    T_range::Tuple{Float64, Float64}  # Temperature range (MeV)
    μ_range::Tuple{Float64, Float64}  # Chemical potential range (MeV)
    n_T_points::Int                   # Number of temperature points
    n_μ_points::Int                   # Number of chemical potential points
    
    # Polyakov potential parameters
    a0::Float64
    a1::Float64
    a2::Float64
    a3::Float64
    b3::Float64
    b4::Float64
    
    function ThermoParams(T_range, μ_range, n_T_points, n_μ_points,
                         a0, a1, a2, a3, b3, b4)
        @assert T_range[1] >= 0 "Minimum temperature must be non-negative"
        @assert T_range[2] > T_range[1] "Temperature range must be valid"
        @assert n_T_points > 0 "Number of T points must be positive"
        @assert n_μ_points > 0 "Number of μ points must be positive"
        new(T_range, μ_range, n_T_points, n_μ_points, a0, a1, a2, a3, b3, b4)
    end
end
```

### Parameter Format Conversion

```julia
"""
Support multiple parameter formats for flexibility.
"""

# Convert NamedTuple to QuarkParams
function QuarkParams(nt::NamedTuple)
    return QuarkParams(
        nt.n_flavors,
        nt.current_masses,
        nt.coupling_G,
        nt.coupling_K,
        nt.cutoff
    )
end

# Convert Dict to QuarkParams
function QuarkParams(d::Dict{Symbol, Any})
    return QuarkParams(
        d[:n_flavors],
        d[:current_masses],
        d[:coupling_G],
        d[:coupling_K],
        d[:cutoff]
    )
end

# Convert QuarkParams to Dict
function Base.Dict(params::QuarkParams)
    return Dict{Symbol, Any}(
        :n_flavors => params.n_flavors,
        :current_masses => params.current_masses,
        :coupling_G => params.coupling_G,
        :coupling_K => params.coupling_K,
        :cutoff => params.cutoff
    )
end
```


### Data Flow

The typical data flow through the system:

```
1. User creates model via factory or direct construction
   ↓
2. Model validates parameters during construction
   ↓
3. User calls high-level calculation function
   ↓
4. Function dispatches to model-specific implementations
   ↓
5. Model implementations use base utilities for:
   - Numerical integration
   - Polyakov loop calculations
   - Safe mathematical operations
   - Distribution functions
   ↓
6. Results returned to user
```

Example calculation flow:

```julia
# 1. Create model
params = QuarkParams(2, [5.5, 5.5], 10.08, 12.36, 650.0)
thermo = ThermoParams((0.0, 300.0), (-500.0, 500.0), 100, 100, 
                      6.75, -1.95, 2.625, -7.44, 0.75, 7.5)
model = create_model(:PNJL, params=params, thermo_params=thermo)

# 2. Calculate thermodynamics at specific point
T = 150.0  # MeV
μ = 300.0  # MeV
result = calculate_thermodynamics(model, T, μ)

# 3. Internal flow:
#    - calculate_thermodynamics calls vacuum_contribution(model)
#    - vacuum_contribution uses vacuum_integral from base utilities
#    - calculate_thermodynamics calls thermal_contribution(model, T, μ)
#    - thermal_contribution solves mass gap using solve_mass_gap
#    - solve_mass_gap uses integrate_momentum from base utilities
#    - Results combined and returned
```


## Correctness Properties

*A property is a characteristic or behavior that should hold true across all valid executions of a system—essentially, a formal statement about what the system should do. Properties serve as the bridge between human-readable specifications and machine-verifiable correctness guarantees.*

The following properties define the correctness criteria for the QCD model architecture refactoring. Each property is universally quantified and references the specific requirements it validates.

### Interface Compliance Properties

**Property 1: All concrete models implement required interface methods**

*For any* concrete model instance that inherits from AbstractQCDModel, calling any of the required interface methods (vacuum_contribution, thermal_contribution, polyakov_potential, dispersion_relation, calculate_mass_vec, calculate_chiral) should not throw a MethodError.

**Validates: Requirements 1.3**


### Base Utilities Properties

**Property 2: Numerical integration produces finite results**

*For any* valid energy function callback and valid temperature/chemical potential parameters, the numerical integration utilities should return a finite Float64 value (not NaN or Inf).

**Validates: Requirements 2.1**

**Property 3: Polyakov log calculations respect physical constraints**

*For any* Polyakov loop values phi and phi_bar in the range [0, 1], the polyakov_log_combination function should return a finite value and handle edge cases (phi=0, phi=1) gracefully.

**Validates: Requirements 2.2**

**Property 4: Safe math functions handle edge cases without crashing**

*For any* input value including edge cases (zero, negative, very large, very small), safe numerical functions (safe_log, safe_vandermonde, safe_divide) should either return a valid result or throw a descriptive exception, but never crash or return NaN unexpectedly.

**Validates: Requirements 2.3, 5.1, 5.2**

**Property 5: Distribution functions return values in valid range**

*For any* valid energy, temperature, chemical potential, and Polyakov loop values, the quark and antiquark distribution functions should return values in the range [0, 1], representing valid probability distributions.

**Validates: Requirements 2.4**


### PNJL Refactoring Properties

**Property 6: Refactored PNJL produces numerically equivalent results**

*For any* valid temperature, chemical potential, and parameter set, the refactored PNJLModel should produce thermodynamic quantities (vacuum contribution, thermal contribution, Polyakov potential, constituent masses, chiral condensates) that are numerically equivalent to the original implementation within floating-point precision (relative error < 1e-12).

**Validates: Requirements 3.2**

**Property 7: Legacy API produces equivalent results through compatibility layer**

*For any* valid legacy API call with old parameter format, the result obtained through the compatibility layer should be numerically equivalent to what the original system would produce (relative error < 1e-12).

**Validates: Requirements 3.5, 6.2, 6.5**


### rPNJL Model Properties

**Property 8: rPNJL mass calculation includes eight-quark terms**

*For any* valid temperature, chemical potential, and rPNJL parameters (g1, g2), the constituent quark mass calculated by RPNJLModel should differ from the standard PNJL mass in a way consistent with the eight-quark interaction term (equation 3.31), specifically: when g1 > 0 or g2 > 0, the rPNJL mass should reflect the repulsive interaction contribution.

**Validates: Requirements 4.3**

**Property 9: rPNJL Polyakov potential includes Vandermonde terms**

*For any* valid temperature and Vandermonde parameters (kappa, T0), the Polyakov potential calculated by RPNJLModel should differ from the standard polynomial potential in a way consistent with the Vandermonde term contribution (equations 3.27-3.29).

**Validates: Requirements 4.4**

**Property 10: Vandermonde term exhibits temperature dependence**

*For any* two different temperatures T1 and T2 (with T1 ≠ T2) and fixed Polyakov loop value phi, the Vandermonde term should produce different values due to the b2(T) temperature dependence, demonstrating that the term is not temperature-independent.

**Validates: Requirements 4.5**


### Parameter Management Properties

**Property 11: Parameter format conversion preserves values**

*For any* valid parameter set, converting from one supported format to another and back should preserve all parameter values (round-trip property): QuarkParams → Dict → QuarkParams should yield equivalent parameters, and similarly for NamedTuple conversions.

**Validates: Requirements 12.2, 12.3, 12.5**

**Property 12: Parameter conversion accepts all supported formats**

*For any* valid parameter data in supported formats (struct, NamedTuple, Dict{Symbol, Any}), the conversion utilities should successfully convert to the target format without throwing exceptions.

**Validates: Requirements 6.3**


### Validation and Error Handling Properties

**Property 13: Invalid parameters trigger validation errors**

*For any* parameter value that violates defined constraints (e.g., negative temperature, Polyakov loop outside [0,1], negative mass), model construction or calculation should throw a ModelParameterError with a descriptive message indicating which constraint was violated.

**Validates: Requirements 14.2, 14.3**

**Property 14: Physical constraints are enforced during calculations**

*For any* calculation that produces intermediate or final values subject to physical constraints (e.g., Polyakov loop in [0,1], non-negative masses, finite energies), the system should validate these constraints and throw PhysicalConstraintError if violated.

**Validates: Requirements 14.4**

**Property 15: All exceptions include descriptive messages**

*For any* custom exception thrown by the system (ModelParameterError, NumericalConvergenceError, PhysicalConstraintError, ModelConfigurationError, NumericalInstabilityError), the exception should contain a non-empty, descriptive message that includes context about what went wrong and suggested remediation.

**Validates: Requirements 13.2**


## Error Handling

### Exception Hierarchy and Responsibilities

The system implements a layered error handling strategy with clear responsibilities at each level:

**Layer 1: Base Utilities**
- Catches numerical library exceptions (integration failures, solver divergence)
- Converts to domain-specific exceptions (NumericalConvergenceError, NumericalInstabilityError)
- Implements fallback strategies (increase quadrature nodes, relax tolerances)
- Logs warnings for recoverable issues

**Layer 2: Model Interface**
- Validates input parameters against constraints
- Throws ModelParameterError for invalid inputs
- Throws PhysicalConstraintError for constraint violations
- Ensures all required methods are implemented

**Layer 3: High-Level API**
- Provides user-friendly error messages
- Catches and re-wraps lower-level exceptions with context
- Logs errors at appropriate levels
- Returns meaningful error codes or exceptions to user scripts

### Fallback Strategies

**Integration Convergence Failures**:
```julia
# Strategy 1: Increase quadrature nodes
try
    result = integrate_momentum(energy_func, T, μ; n_nodes=64)
catch e
    @warn "Integration failed, doubling nodes" exception=e
    result = integrate_momentum(energy_func, T, μ; n_nodes=128)
end
```

**Solver Convergence Failures**:
```julia
# Strategy 1: Relax tolerance
# Strategy 2: Adjust initial guess
# Strategy 3: Try different solver algorithm
# Maximum 3 retry attempts before throwing NumericalConvergenceError
```

**Negative Vandermonde Discriminant**:
```julia
# Log warning but continue calculation
# Physical interpretation may still be valid
@warn "Negative Vandermonde discriminant" phi=phi b2=b2 discriminant=discriminant
return calculate_vandermonde_term(phi, b2, discriminant)
```

### Error Message Guidelines

All error messages should follow this structure:
1. **What went wrong**: Clear description of the error
2. **Context**: Relevant parameter values and state
3. **Why it matters**: Impact of the error
4. **What to do**: Suggested remediation steps

Example:
```julia
throw(ModelParameterError(
    "temperature",
    T,
    "0.0 <= T <= 500.0",
    "Temperature $T MeV is outside physical range. " *
    "QCD phase transitions occur in 0-300 MeV range. " *
    "Check input parameters or use extended range if studying extreme conditions."
))
```


## Testing Strategy

### Dual Testing Approach

The testing strategy employs both unit tests and property-based tests as complementary approaches:

**Unit Tests**: Focus on specific examples, edge cases, and integration points
- Specific parameter combinations with known results
- Edge cases (T=0, μ=0, extreme values)
- Error conditions and exception handling
- Integration between components
- Regression tests comparing to original implementation

**Property-Based Tests**: Verify universal properties across randomized inputs
- Interface compliance across all model types
- Numerical stability across parameter ranges
- Physical constraint satisfaction
- Mathematical properties (distribution ranges, convergence)
- Round-trip properties (parameter conversion)

Together, these approaches provide comprehensive coverage: unit tests catch concrete bugs in specific scenarios, while property tests verify general correctness across the input space.

### Property-Based Testing Configuration

**Library**: Use Julia's `PropertyBasedTesting.jl` or `Supposition.jl` for property-based testing

**Configuration**:
- Minimum 100 iterations per property test (due to randomization)
- Each test tagged with feature name and property number
- Tag format: `# Feature: qcd-model-architecture-refactoring, Property N: [property description]`

**Example Property Test**:
```julia
using PropertyBasedTesting

@testset "Property 6: PNJL numerical equivalence" begin
    # Feature: qcd-model-architecture-refactoring, Property 6: Refactored PNJL produces numerically equivalent results
    @check iterations=100 function pnjl_numerical_equivalence(
        T = Float64(0.0..300.0),
        μ = Float64(-500.0..500.0)
    )
        # Create refactored model
        params = get_test_pnjl_params()
        model_new = PNJLModel(params.quark, params.thermo, 650.0)
        
        # Calculate with new implementation
        result_new = calculate_thermodynamics(model_new, T, μ)
        
        # Calculate with original implementation
        result_old = calculate_thermodynamics_original(T, μ, params)
        
        # Check numerical equivalence
        relative_error = abs(result_new - result_old) / (abs(result_old) + 1e-10)
        return relative_error < 1e-12
    end
end
```


### Test Organization

```
tests/
├── unit/
│   ├── models/
│   │   ├── test_pnjl.jl              # PNJL model unit tests
│   │   ├── test_rpnjl.jl             # rPNJL model unit tests
│   │   ├── test_exceptions.jl        # Exception handling tests
│   │   └── test_validation.jl        # Parameter validation tests
│   ├── base/
│   │   ├── test_integrals.jl         # Integration utilities tests
│   │   ├── test_polyakov.jl          # Polyakov calculations tests
│   │   ├── test_safe_math.jl         # Safe math operations tests
│   │   ├── test_distributions.jl     # Distribution functions tests
│   │   └── test_mass_chiral.jl       # Mass/chiral calculations tests
│   └── compatibility/
│       └── test_legacy_api.jl        # Backward compatibility tests
├── property/
│   ├── test_interface_compliance.jl  # Property 1
│   ├── test_base_utilities.jl        # Properties 2-5
│   ├── test_pnjl_equivalence.jl      # Properties 6-7
│   ├── test_rpnjl_physics.jl         # Properties 8-10
│   ├── test_parameters.jl            # Properties 11-12
│   └── test_validation.jl            # Properties 13-15
├── regression/
│   └── pnjl_refactoring/
│       ├── test_thermodynamics.jl    # Compare all thermodynamic quantities
│       ├── test_phase_diagram.jl     # Compare phase transition points
│       └── baseline_data/            # Reference data from original implementation
├── integration/
│   ├── test_full_workflow.jl         # End-to-end calculation workflows
│   └── test_model_switching.jl       # Switching between models
├── perf/
│   └── models/
│       ├── benchmark_pnjl.jl         # PNJL performance benchmarks
│       └── benchmark_rpnjl.jl        # rPNJL performance benchmarks
└── edge_cases/
    ├── test_zero_temperature.jl      # T=0 edge cases
    ├── test_extreme_mu.jl            # Extreme chemical potential
    └── test_numerical_limits.jl      # Numerical boundary conditions
```

### Test Coverage Requirements

- **Core functionality**: ≥80% code coverage
- **Base utilities**: 100% coverage (critical numerical code)
- **Error handling**: All exception types must be triggered in tests
- **Property tests**: All 15 correctness properties must have corresponding tests
- **Regression tests**: All thermodynamic quantities compared against baseline

### Continuous Integration

Tests should be run in CI pipeline with:
1. Unit tests on every commit
2. Property tests on every pull request
3. Regression tests before merging to main
4. Performance benchmarks weekly
5. Coverage reports generated and tracked


## Implementation Approach

### Phased Implementation Strategy

The refactoring will be implemented in phases to minimize risk and enable incremental validation:

**Phase 1: Foundation (Week 1)**
- Create core interface definitions (AbstractQCDModel, domain abstracts)
- Implement exception hierarchy
- Set up base utilities (integrals, safe math, distributions)
- Establish testing infrastructure
- **Validation**: Base utilities pass unit tests

**Phase 2: PNJL Refactoring (Week 1-2)**
- Extract PNJL implementation into new architecture
- Implement PNJLModel with all required methods
- Create regression test suite comparing to original
- Implement compatibility layer for legacy API
- **Validation**: All regression tests pass with <1e-12 error

**Phase 3: rPNJL Implementation (Week 2-3)**
- Implement RPNJLModel with eight-quark terms
- Implement Vandermonde-modified Polyakov potential
- Add rPNJL-specific parameter validation
- Create rPNJL test suite
- **Validation**: rPNJL physics tests pass, eight-quark and Vandermonde terms verified

**Phase 4: Integration and Optimization (Week 3)**
- Implement factory pattern
- Add parameter conversion utilities
- Performance optimization and profiling
- Complete documentation
- **Validation**: Performance within 5% of baseline, all tests pass

### Development Workflow

1. **Feature Branch**: Create branch for each phase
2. **Test-Driven**: Write tests before implementation
3. **Incremental Commits**: Small, focused commits with clear messages
4. **Code Review**: All changes reviewed before merging
5. **Continuous Testing**: Run test suite on every commit

### Migration Path for Users

Users can migrate to the new architecture gradually:

**Option 1: Immediate Migration**
```julia
# Old code
result = pnjl_thermodynamics(T, μ)

# New code
model = create_model(:PNJL, params=default_params)
result = calculate_thermodynamics(model, T, μ)
```

**Option 2: Compatibility Layer (No Changes Required)**
```julia
# Existing code continues to work
result = pnjl_thermodynamics(T, μ)  # Uses compatibility layer internally
```

**Option 3: Gradual Migration**
```julia
# Mix old and new APIs during transition
model = create_model(:PNJL, params=convert_legacy_params(old_params))
result = calculate_thermodynamics(model, T, μ)
```


## Risk Analysis and Mitigation

### Technical Risks

**Risk 1: Numerical Precision Loss**
- **Probability**: Medium
- **Impact**: High (invalidates scientific results)
- **Mitigation**: 
  - Comprehensive regression testing with tight tolerances (<1e-12)
  - Bit-for-bit comparison where possible
  - Use same numerical libraries and algorithms as original
  - Property tests to verify mathematical properties
- **Contingency**: If precision loss detected, investigate and fix before proceeding

**Risk 2: Performance Degradation**
- **Probability**: Medium
- **Impact**: Medium (slower simulations)
- **Mitigation**:
  - Performance benchmarks at each phase
  - Use @inline hints for critical paths
  - Leverage Julia's type specialization
  - Profile and optimize hot paths
  - Target: <5% degradation
- **Contingency**: If >5% degradation, profile and optimize before proceeding

**Risk 3: Breaking Changes to User Scripts**
- **Probability**: Low (with compatibility layer)
- **Impact**: High (disrupts research workflows)
- **Mitigation**:
  - Comprehensive compatibility layer
  - Test with real user scripts
  - Maintain all public exports
  - Clear migration documentation
- **Contingency**: Extend compatibility layer to cover additional use cases

**Risk 4: Incomplete Interface Implementation**
- **Probability**: Low
- **Impact**: Medium (runtime errors)
- **Mitigation**:
  - Property tests for interface compliance
  - Clear error messages for missing methods
  - Comprehensive documentation of required methods
  - Example implementations for reference
- **Contingency**: Add runtime checks and better error messages

**Risk 5: rPNJL Physics Implementation Errors**
- **Probability**: Medium (new model, complex equations)
- **Impact**: High (incorrect physics)
- **Mitigation**:
  - Careful implementation of equations 3.27-3.31
  - Unit tests for each term separately
  - Property tests for expected physical behavior
  - Comparison with published results if available
  - Code review by physics expert
- **Contingency**: Extensive validation against literature and known limits

### Project Risks

**Risk 6: Timeline Overrun**
- **Probability**: Medium
- **Impact**: Medium
- **Mitigation**:
  - Phased approach allows early delivery of core functionality
  - Clear milestones and validation criteria
  - Regular progress reviews
- **Contingency**: Prioritize P0 requirements, defer optional features

**Risk 7: Insufficient Testing**
- **Probability**: Low
- **Impact**: High (bugs in production)
- **Mitigation**:
  - Dual testing strategy (unit + property tests)
  - 80% coverage requirement
  - Regression test suite
  - CI/CD pipeline
- **Contingency**: Extend testing phase if coverage insufficient

### Success Criteria

The refactoring is considered successful when:
1. ✅ All regression tests pass with <1e-12 relative error
2. ✅ Performance degradation <5% on critical paths
3. ✅ All 15 correctness properties verified by property tests
4. ✅ Test coverage ≥80% for core functionality
5. ✅ All existing user scripts work without modification (via compatibility layer)
6. ✅ rPNJL model produces physically reasonable results
7. ✅ Documentation complete and reviewed
8. ✅ Code review completed for all components


## Logging and Monitoring

### Logging Strategy

The system uses Julia's standard `Logging` library with four levels:

**DEBUG Level**: Detailed internal state for debugging
```julia
@debug "Mass gap iteration" iteration=iter M_old=M M_new=M_new residual=residual
@debug "Integration parameters" n_nodes=n_nodes cutoff=cutoff T=T μ=μ
```

**INFO Level**: Normal operation milestones
```julia
@info "Starting thermodynamic scan" T_range=T_range μ_range=μ_range n_points=n_points
@info "Model created" model_type=typeof(model) parameters=summary(params)
@info "Calculation completed" elapsed_time=elapsed convergence=converged
```

**WARN Level**: Recoverable issues and performance concerns
```julia
@warn "Integration convergence slow" iterations=iter tolerance=tol
@warn "Negative Vandermonde discriminant" phi=phi b2=b2 discriminant=discriminant
@warn "Performance degradation detected" expected_time=t_expected actual_time=t_actual
```

**ERROR Level**: Serious errors (program continues)
```julia
@error "Parameter validation failed" parameter=param value=value constraint=constraint
@error "Solver diverged" max_iterations=max_iter residual=residual
```

### Structured Logging Format

Use key-value pairs for machine-readable logs:
```julia
@info "Operation completed" operation="thermodynamic_scan" status="success" 
      duration_ms=elapsed_ms n_points=n_points convergence_rate=0.95
```

### Performance Monitoring

Key metrics to monitor:
- **Calculation time**: Track time for each major operation
- **Convergence iterations**: Monitor solver convergence behavior
- **Memory allocation**: Track allocations in hot paths
- **Cache hit rates**: Monitor effectiveness of memoization

Example monitoring:
```julia
function calculate_thermodynamics(model::AbstractQCDModel, T::Float64, μ::Float64)
    start_time = time()
    
    result = # ... calculation ...
    
    elapsed = time() - start_time
    @info "Thermodynamics calculated" model=typeof(model) T=T μ=μ elapsed_ms=elapsed*1000
    
    return result
end
```

### Log Configuration

Users can configure logging:
```julia
# Set global log level
using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

# Set per-module log level
ENV["JULIA_DEBUG"] = "QCDModels"  # Enable DEBUG for QCDModels module

# Log to file
using LoggingExtras
file_logger = FileLogger("qcd_calculations.log")
global_logger(file_logger)
```

### Monitoring Best Practices

1. **Avoid logging large arrays**: Use summary statistics instead
   ```julia
   @debug "Mass vector calculated" n_flavors=length(M_vec) mean_mass=mean(M_vec) max_mass=maximum(M_vec)
   ```

2. **Use structured data**: Key-value pairs for easy parsing
3. **Log at appropriate levels**: Don't spam INFO with DEBUG information
4. **Include context**: Always include relevant parameters in log messages
5. **Performance-aware**: Logging should not significantly impact performance


## Summary

This design document specifies a comprehensive refactoring of the QCD model library from a monolithic PNJL implementation to a polymorphic, extensible architecture. The key design decisions are:

1. **Multi-dispatch polymorphism**: Leverages Julia's type system for zero-overhead extensibility
2. **Layered architecture**: Clear separation between interfaces, implementations, and utilities
3. **Backward compatibility**: Existing scripts continue working through compatibility layer
4. **Safety-first**: Comprehensive validation, error handling, and safe numerical operations
5. **Test-driven**: Dual testing strategy with unit tests and property-based tests
6. **Phased implementation**: Incremental delivery with validation at each phase

The architecture supports the immediate implementation of PNJL and rPNJL models while providing a clear path for adding future models (magnetic PNJL, etc.) without modifying core infrastructure. The design maintains numerical precision, meets performance requirements, and provides comprehensive error handling and logging for production use.

