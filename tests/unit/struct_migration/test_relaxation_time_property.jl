"""
Property-based tests for RelaxationTime module struct equivalence.

Tests Property 1: Struct-NamedTuple Equivalence (RelaxationTime)
Validates: Requirements 1.2, 11.1

This test verifies that relaxation_times produces identical results
whether called with QuarkParams/ThermoParams structs or NamedTuple parameters.
"""

using Test
using Supposition

# Load required modules
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src/relaxtime"))

include("../../../src/Constants_PNJL.jl")
include("../../../src/relaxtime/RelaxationTime.jl")
include("../../../src/relaxtime/EffectiveCouplings.jl")
include("../../../src/relaxtime/OneLoopIntegrals.jl")
include("../../../src/integration/GaussLegendre.jl")

using .RelaxationTime
using .Constants_PNJL
using .EffectiveCouplings
using .OneLoopIntegrals: A
using .GaussLegendre: gauleg

# Load test utilities
include("test_utils.jl")
using .Main: QuarkParams, ThermoParams, as_namedtuple, approx_equal

@testset "RelaxationTime Property Tests" begin
    
    # ========================================================================
    # Property 1: Struct-NamedTuple Equivalence for relaxation_times
    # ========================================================================
    # **Validates: Requirements 1.2, 11.1**
    #
    # This property test verifies that relaxation_times produces
    # numerically identical results when called with struct parameters vs
    # NamedTuple parameters across a wide range of random inputs.
    
    @testset "Property: Struct-NamedTuple Equivalence" begin
        println("\n" * "="^70)
        println("Property Test: Relaxation Times Struct-NamedTuple Equivalence")
        println("Testing relaxation_times with random parameters")
        println("Note: Using pre-computed rates to avoid expensive computation")
        println("="^70)
        
        RTOL = 1e-12  # Relative tolerance for floating-point comparison
        ATOL = 1e-14  # Absolute tolerance
        
        # Pre-compute Gauss-Legendre nodes and weights for A function
        nodes_p, weights_p = gauleg(0.0, 20.0, 64)
        
        @check max_examples=5 function property_relaxation_times_equivalence(
            # Generate random quark masses (in fm⁻¹)
            m_u = Data.Floats{Float64}(minimum=0.5, maximum=2.0),
            m_s = Data.Floats{Float64}(minimum=2.0, maximum=5.0),
            # Generate random chemical potentials (in fm⁻¹)
            μ_u = Data.Floats{Float64}(minimum=0.0, maximum=0.5),
            μ_s = Data.Floats{Float64}(minimum=0.0, maximum=0.5),
            # Generate random thermodynamic parameters
            T = Data.Floats{Float64}(minimum=0.1, maximum=0.3),
            Φ = Data.Floats{Float64}(minimum=0.1, maximum=0.9),
            Φbar = Data.Floats{Float64}(minimum=0.1, maximum=0.9),
            # Generate random densities
            n_u = Data.Floats{Float64}(minimum=0.01, maximum=0.5),
            n_s = Data.Floats{Float64}(minimum=0.01, maximum=0.3),
        )
            # Compute A functions for both u and s quarks
            A_u = A(m_u, μ_u, T, Φ, Φbar, nodes_p, weights_p)
            A_s = A(m_s, μ_s, T, Φ, Φbar, nodes_p, weights_p)
            
            # Compute G functions
            G_u = calculate_G_from_A(A_u, m_u)
            G_s = calculate_G_from_A(A_s, m_s)
            
            # Compute K_coeffs
            G_fm2 = Constants_PNJL.G_fm2
            K_fm5 = Constants_PNJL.K_fm5
            K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
            
            # Create struct parameters
            q_struct = QuarkParams(
                (u=m_u, d=m_u, s=m_s),
                (u=μ_u, d=μ_u, s=μ_s)
            )
            t_struct = ThermoParams(T, Φ, Φbar, 0.0)
            
            # Create NamedTuple parameters
            q_nt = (
                m = (u=m_u, d=m_u, s=m_s),
                μ = (u=μ_u, d=μ_u, s=μ_s)
            )
            t_nt = (T=T, Φ=Φ, Φbar=Φbar, ξ=0.0)
            
            # Create densities (same for both tests)
            densities = (
                u = n_u,
                d = n_u,  # isospin symmetry
                s = n_s,
                ubar = n_u,
                dbar = n_u,
                sbar = n_s
            )
            
            # Pre-compute rates once to avoid expensive computation in the test
            # This tests that the struct interface works correctly with pre-computed rates
            rates = compute_average_rates(
                q_nt, t_nt, K_coeffs,
                p_nodes=2,
                angle_nodes=2,
                phi_nodes=2,
                n_sigma_points=3
            )
            
            # Compute relaxation times with struct parameters using pre-computed rates
            result_struct = relaxation_times(
                q_struct, t_struct, K_coeffs;
                densities=densities,
                existing_rates=rates
            )
            
            # Compute relaxation times with NamedTuple parameters using pre-computed rates
            result_nt = relaxation_times(
                q_nt, t_nt, K_coeffs;
                densities=densities,
                existing_rates=rates
            )
            
            # Verify tau values are equivalent
            for flavor in [:u, :d, :s, :ubar, :dbar, :sbar]
                tau_struct = getproperty(result_struct.tau, flavor)
                tau_nt = getproperty(result_nt.tau, flavor)
                
                # Both should be positive and finite
                tau_struct > 0.0 || error("tau.$flavor must be positive (struct): $tau_struct")
                tau_nt > 0.0 || error("tau.$flavor must be positive (nt): $tau_nt")
                isfinite(tau_struct) || error("tau.$flavor must be finite (struct): $tau_struct")
                isfinite(tau_nt) || error("tau.$flavor must be finite (nt): $tau_nt")
                
                # Verify numerical equivalence
                isapprox(tau_struct, tau_nt, rtol=RTOL, atol=ATOL) ||
                    error("tau.$flavor mismatch: struct=$tau_struct, nt=$tau_nt, " *
                          "diff=$(abs(tau_struct - tau_nt))")
            end
            
            # Verify tau_inv values are equivalent
            for flavor in [:u, :d, :s, :ubar, :dbar, :sbar]
                tau_inv_struct = getproperty(result_struct.tau_inv, flavor)
                tau_inv_nt = getproperty(result_nt.tau_inv, flavor)
                
                # Both should be non-negative and finite
                tau_inv_struct >= 0.0 || error("tau_inv.$flavor must be non-negative (struct): $tau_inv_struct")
                tau_inv_nt >= 0.0 || error("tau_inv.$flavor must be non-negative (nt): $tau_inv_nt")
                isfinite(tau_inv_struct) || error("tau_inv.$flavor must be finite (struct): $tau_inv_struct")
                isfinite(tau_inv_nt) || error("tau_inv.$flavor must be finite (nt): $tau_inv_nt")
                
                # Verify numerical equivalence
                isapprox(tau_inv_struct, tau_inv_nt, rtol=RTOL, atol=ATOL) ||
                    error("tau_inv.$flavor mismatch: struct=$tau_inv_struct, nt=$tau_inv_nt")
            end
            
            # Return true to indicate test passed
            true
        end
        
        println("✓ Property test passed: Struct-NamedTuple equivalence verified for relaxation_times")
    end
    
    # ========================================================================
    # Property 2: Struct-NamedTuple Equivalence for compute_average_rates
    # ========================================================================
    
    @testset "Property: compute_average_rates Equivalence" begin
        println("\n" * "="^70)
        println("Property Test: compute_average_rates Struct-NamedTuple Equivalence")
        println("Testing compute_average_rates with random parameters")
        println("Note: Testing with minimal integration nodes for performance")
        println("="^70)
        
        RTOL = 1e-12
        ATOL = 1e-14
        
        nodes_p, weights_p = gauleg(0.0, 20.0, 64)
        
        @check max_examples=3 function property_compute_average_rates_equivalence(
            m_u = Data.Floats{Float64}(minimum=0.5, maximum=2.0),
            m_s = Data.Floats{Float64}(minimum=2.0, maximum=5.0),
            μ_u = Data.Floats{Float64}(minimum=0.0, maximum=0.5),
            μ_s = Data.Floats{Float64}(minimum=0.0, maximum=0.5),
            T = Data.Floats{Float64}(minimum=0.1, maximum=0.3),
            Φ = Data.Floats{Float64}(minimum=0.1, maximum=0.9),
            Φbar = Data.Floats{Float64}(minimum=0.1, maximum=0.9),
        )
            # Setup parameters
            A_u = A(m_u, μ_u, T, Φ, Φbar, nodes_p, weights_p)
            A_s = A(m_s, μ_s, T, Φ, Φbar, nodes_p, weights_p)
            G_u = calculate_G_from_A(A_u, m_u)
            G_s = calculate_G_from_A(A_s, m_s)
            K_coeffs = calculate_effective_couplings(
                Constants_PNJL.G_fm2, Constants_PNJL.K_fm5, G_u, G_s
            )
            
            q_struct = QuarkParams((u=m_u, d=m_u, s=m_s), (u=μ_u, d=μ_u, s=μ_s))
            t_struct = ThermoParams(T, Φ, Φbar, 0.0)
            q_nt = (m=(u=m_u, d=m_u, s=m_s), μ=(u=μ_u, d=μ_u, s=μ_s))
            t_nt = (T=T, Φ=Φ, Φbar=Φbar, ξ=0.0)
            
            # Compute average rates with struct parameters (minimal nodes)
            rates_struct = compute_average_rates(
                q_struct, t_struct, K_coeffs,
                p_nodes=2,
                angle_nodes=2,
                phi_nodes=2,
                n_sigma_points=3
            )
            
            # Compute average rates with NamedTuple parameters
            rates_nt = compute_average_rates(
                q_nt, t_nt, K_coeffs,
                p_nodes=2,
                angle_nodes=2,
                phi_nodes=2,
                n_sigma_points=3
            )
            
            # Verify all rates are equivalent (test a subset for performance)
            test_processes = [:uu_to_uu, :ss_to_ss, :ud_to_ud, :us_to_us]
            for process in test_processes
                if hasproperty(rates_struct, process) && hasproperty(rates_nt, process)
                    rate_struct = getproperty(rates_struct, process)
                    rate_nt = getproperty(rates_nt, process)
                    
                    # Both should be non-negative and finite
                    rate_struct >= 0.0 || error("rates.$process must be non-negative (struct): $rate_struct")
                    rate_nt >= 0.0 || error("rates.$process must be non-negative (nt): $rate_nt")
                    isfinite(rate_struct) || error("rates.$process must be finite (struct): $rate_struct")
                    isfinite(rate_nt) || error("rates.$process must be finite (nt): $rate_nt")
                    
                    # Verify numerical equivalence
                    isapprox(rate_struct, rate_nt, rtol=RTOL, atol=ATOL) ||
                        error("rates.$process mismatch: struct=$rate_struct, nt=$rate_nt")
                end
            end
            
            # Return true to indicate test passed
            true
        end
        
        println("✓ Property test passed: compute_average_rates equivalence verified")
    end
    
end

println("\n" * "="^70)
println("RelaxationTime Property Tests Complete!")
println("="^70)
println("\nAll property tests passed:")
println("  1. Struct-NamedTuple equivalence for relaxation_times")
println("  2. Struct-NamedTuple equivalence for compute_average_rates")
println("  3. Physical constraints maintained (positive tau, non-negative rates)")
println("  4. All values are finite")
println("="^70)
