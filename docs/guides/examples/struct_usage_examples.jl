# Parameter Struct Usage Examples
#
# This file demonstrates how to use QuarkParams and ThermoParams structs
# in the PNJL model codebase. These examples show the recommended patterns
# for new code.

using Main.ParameterTypes: QuarkParams, ThermoParams, as_namedtuple

# Include necessary modules (adjust paths as needed)
include("../../src/relaxtime/RelaxationTime.jl")
using .RelaxationTime

# =============================================================================
# Example 1: Basic Struct Usage (Recommended Pattern)
# =============================================================================

function example_1_basic_struct_usage()
    println("\n=== Example 1: Basic Struct Usage ===\n")
    
    # Create QuarkParams struct with type safety
    quark_params = QuarkParams(
        m=(u=1.52, d=1.52, s=3.04),  # Quark masses in fm⁻¹
        μ=(u=0.3, d=0.3, s=0.3)      # Chemical potentials in fm⁻¹
    )
    
    # Create ThermoParams struct
    thermo_params = ThermoParams(
        0.15,  # Temperature T in fm⁻¹
        0.5,   # Polyakov loop Φ
        0.5,   # Conjugate Polyakov loop Φbar
        0.0    # Anisotropy parameter ξ
    )
    
    # Define coupling coefficients (still uses NamedTuple)
    K_coeffs = (
        K_pi=3.67,
        K_K=5.50,
        K_eta=0.0,
        K_sigma=3.67,
        K_zeta=0.0,
        K_delta=0.0,
        K_kappa=5.50
    )
    
    # Define quark densities
    densities = (
        u=0.1,
        d=0.1,
        s=0.05,
        ubar=0.08,
        dbar=0.08,
        sbar=0.04
    )
    
    # Calculate relaxation times using struct parameters
    result = relaxation_times(
        quark_params,
        thermo_params,
        K_coeffs;
        densities=densities
    )
    
    println("Relaxation times (using structs):")
    println("  τ_u = ", result.tau.u, " fm")
    println("  τ_s = ", result.tau.s, " fm")
    
    return result
end

# =============================================================================
# Example 2: NamedTuple Usage (Backward Compatible)
# =============================================================================

function example_2_namedtuple_usage()
    println("\n=== Example 2: NamedTuple Usage (Backward Compatible) ===\n")
    
    # Create parameters using NamedTuples (old style, still works!)
    quark_params = (
        m=(u=1.52, d=1.52, s=3.04),
        μ=(u=0.3, d=0.3, s=0.3)
    )
    
    thermo_params = (
        T=0.15,
        Φ=0.5,
        Φbar=0.5,
        ξ=0.0
    )
    
    K_coeffs = (
        K_pi=3.67,
        K_K=5.50,
        K_eta=0.0,
        K_sigma=3.67,
        K_zeta=0.0,
        K_delta=0.0,
        K_kappa=5.50
    )
    
    densities = (
        u=0.1,
        d=0.1,
        s=0.05,
        ubar=0.08,
        dbar=0.08,
        sbar=0.04
    )
    
    # Calculate relaxation times using NamedTuple parameters
    # This produces IDENTICAL results to Example 1
    result = relaxation_times(
        quark_params,
        thermo_params,
        K_coeffs;
        densities=densities
    )
    
    println("Relaxation times (using NamedTuples):")
    println("  τ_u = ", result.tau.u, " fm")
    println("  τ_s = ", result.tau.s, " fm")
    
    return result
end

# =============================================================================
# Example 3: Mixed Usage (Struct + NamedTuple)
# =============================================================================

function example_3_mixed_usage()
    println("\n=== Example 3: Mixed Usage ===\n")
    
    # Use struct for quark parameters (type safety)
    quark_params = QuarkParams(
        m=(u=1.52, d=1.52, s=3.04),
        μ=(u=0.3, d=0.3, s=0.3)
    )
    
    # Use NamedTuple for thermo parameters (quick prototyping)
    thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
    
    K_coeffs = (
        K_pi=3.67,
        K_K=5.50,
        K_eta=0.0,
        K_sigma=3.67,
        K_zeta=0.0,
        K_delta=0.0,
        K_kappa=5.50
    )
    
    densities = (
        u=0.1,
        d=0.1,
        s=0.05,
        ubar=0.08,
        dbar=0.08,
        sbar=0.04
    )
    
    # Mix struct and NamedTuple - both work together!
    result = relaxation_times(
        quark_params,      # Struct
        thermo_params,     # NamedTuple
        K_coeffs;
        densities=densities
    )
    
    println("Relaxation times (mixed struct/NamedTuple):")
    println("  τ_u = ", result.tau.u, " fm")
    println("  τ_s = ", result.tau.s, " fm")
    
    return result
end

# =============================================================================
# Example 4: Converting Between Formats
# =============================================================================

function example_4_conversion()
    println("\n=== Example 4: Converting Between Formats ===\n")
    
    # Start with a struct
    quark_struct = QuarkParams(
        m=(u=1.52, d=1.52, s=3.04),
        μ=(u=0.3, d=0.3, s=0.3)
    )
    
    println("Original struct:")
    println("  ", quark_struct)
    
    # Convert struct to NamedTuple
    quark_nt = as_namedtuple(quark_struct)
    
    println("\nConverted to NamedTuple:")
    println("  ", quark_nt)
    
    # Convert back to struct
    quark_struct_2 = QuarkParams(quark_nt)
    
    println("\nConverted back to struct:")
    println("  ", quark_struct_2)
    
    # Verify they're equivalent
    println("\nAre they equivalent? ", quark_struct == quark_struct_2)
    
    return quark_struct, quark_nt, quark_struct_2
end

# =============================================================================
# Example 5: Migrating Existing Code
# =============================================================================

# BEFORE: Old code using NamedTuples
function calculate_transport_old_style(T, μ_B)
    quark_params = (
        m=(u=1.52, d=1.52, s=3.04),
        μ=(u=μ_B/3, d=μ_B/3, s=μ_B/3)
    )
    
    thermo_params = (T=T, Φ=0.5, Φbar=0.5, ξ=0.0)
    
    K_coeffs = (
        K_pi=3.67, K_K=5.50, K_eta=0.0,
        K_sigma=3.67, K_zeta=0.0, K_delta=0.0, K_kappa=5.50
    )
    
    densities = (u=0.1, d=0.1, s=0.05, ubar=0.08, dbar=0.08, sbar=0.04)
    
    return relaxation_times(quark_params, thermo_params, K_coeffs; densities=densities)
end

# AFTER: New code using structs
function calculate_transport_new_style(T, μ_B)
    # Only change: use struct constructors instead of NamedTuple literals
    quark_params = QuarkParams(
        m=(u=1.52, d=1.52, s=3.04),
        μ=(u=μ_B/3, d=μ_B/3, s=μ_B/3)
    )
    
    thermo_params = ThermoParams(T, 0.5, 0.5, 0.0)
    
    K_coeffs = (
        K_pi=3.67, K_K=5.50, K_eta=0.0,
        K_sigma=3.67, K_zeta=0.0, K_delta=0.0, K_kappa=5.50
    )
    
    densities = (u=0.1, d=0.1, s=0.05, ubar=0.08, dbar=0.08, sbar=0.04)
    
    # Function call is IDENTICAL
    return relaxation_times(quark_params, thermo_params, K_coeffs; densities=densities)
end

function example_5_migration()
    println("\n=== Example 5: Migrating Existing Code ===\n")
    
    T = 0.15
    μ_B = 0.8
    
    println("Running old-style function (NamedTuples)...")
    result_old = calculate_transport_old_style(T, μ_B)
    println("  τ_u = ", result_old.tau.u, " fm")
    
    println("\nRunning new-style function (structs)...")
    result_new = calculate_transport_new_style(T, μ_B)
    println("  τ_u = ", result_new.tau.u, " fm")
    
    println("\nResults are identical? ", result_old.tau.u ≈ result_new.tau.u)
    
    return result_old, result_new
end

# =============================================================================
# Example 6: Error Handling with Structs
# =============================================================================

function example_6_error_handling()
    println("\n=== Example 6: Error Handling with Structs ===\n")
    
    # Structs catch errors at construction time
    println("Attempting to create QuarkParams with missing field...")
    try
        # This will fail because μ is missing
        bad_params = QuarkParams((m=(u=1.52, d=1.52, s=3.04),))
    catch e
        println("  Caught error: ", e)
    end
    
    # NamedTuples catch errors at usage time
    println("\nAttempting to use NamedTuple with missing field...")
    try
        bad_params = (m=(u=1.52, d=1.52, s=3.04),)  # Missing μ
        thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
        K_coeffs = (K_pi=3.67, K_K=5.50, K_eta=0.0, K_sigma=3.67, K_zeta=0.0, K_delta=0.0, K_kappa=5.50)
        densities = (u=0.1, d=0.1, s=0.05, ubar=0.08, dbar=0.08, sbar=0.04)
        
        # Error happens here when trying to access μ
        result = relaxation_times(bad_params, thermo_params, K_coeffs; densities=densities)
    catch e
        println("  Caught error: ", e)
    end
    
    println("\nConclusion: Structs catch errors earlier (at construction time)")
end

# =============================================================================
# Example 7: Working with Lower-Level Functions
# =============================================================================

function example_7_lower_level_functions()
    println("\n=== Example 7: Lower-Level Functions ===\n")
    
    using .RelaxationTime.AverageScatteringRate
    
    # Create parameters
    quark_params = QuarkParams(
        m=(u=1.52, d=1.52, s=3.04),
        μ=(u=0.3, d=0.3, s=0.3)
    )
    
    thermo_params = ThermoParams(0.15, 0.5, 0.5, 0.0)
    
    K_coeffs = (
        K_pi=3.67, K_K=5.50, K_eta=0.0,
        K_sigma=3.67, K_zeta=0.0, K_delta=0.0, K_kappa=5.50
    )
    
    # Call lower-level function directly
    println("Computing average scattering rate for u+u → u+u...")
    rate = average_scattering_rate(
        :uu_to_uu,
        quark_params,
        thermo_params,
        K_coeffs
    )
    
    println("  <σv> = ", rate, " fm²")
    
    return rate
end

# =============================================================================
# Main: Run All Examples
# =============================================================================

function run_all_examples()
    println("=" ^ 70)
    println("Parameter Struct Usage Examples")
    println("=" ^ 70)
    
    example_1_basic_struct_usage()
    example_2_namedtuple_usage()
    example_3_mixed_usage()
    example_4_conversion()
    example_5_migration()
    example_6_error_handling()
    example_7_lower_level_functions()
    
    println("\n" * "=" ^ 70)
    println("All examples completed!")
    println("=" ^ 70)
end

# Uncomment to run all examples:
# run_all_examples()
