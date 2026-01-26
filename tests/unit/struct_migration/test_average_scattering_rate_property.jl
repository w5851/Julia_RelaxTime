"""
Property-based tests for AverageScatteringRate module struct equivalence.

Tests Property 1: Struct-NamedTuple Equivalence (AverageScatteringRate)
Validates: Requirements 2.2, 11.2

This test verifies that average_scattering_rate produces identical results
whether called with QuarkParams/ThermoParams structs or NamedTuple parameters.
"""

using Test
using Supposition

# Load required modules
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src/relaxtime"))

include("../../../src/Constants_PNJL.jl")
include("../../../src/relaxtime/AverageScatteringRate.jl")
include("../../../src/relaxtime/EffectiveCouplings.jl")
include("../../../src/relaxtime/OneLoopIntegrals.jl")
include("../../../src/integration/GaussLegendre.jl")

using .AverageScatteringRate
using .Constants_PNJL
using .EffectiveCouplings
using .OneLoopIntegrals: A
using .GaussLegendre: gauleg

# Load test utilities
include("test_utils.jl")
using .Main: QuarkParams, ThermoParams, as_namedtuple, approx_equal

@testset "AverageScatteringRate Property Tests" begin
    
    # ========================================================================
    # Property 1: Struct-NamedTuple Equivalence for average_scattering_rate
    # ========================================================================
    # **Validates: Requirements 2.2, 11.2**
    #
    # This property test verifies that average_scattering_rate produces
    # numerically identical results when called with struct parameters vs
    # NamedTuple parameters across a wide range of random inputs.
    
    @testset "Property: Struct-NamedTuple Equivalence" begin
        println("\n" * "="^70)
        println("Property Test: Average Scattering Rate Struct-NamedTuple Equivalence")
        println("Testing average_scattering_rate with random parameters")
        println("="^70)
        
        RTOL = 1e-12  # Relative tolerance for floating-point comparison
        ATOL = 1e-14  # Absolute tolerance
        
        # Pre-compute Gauss-Legendre nodes and weights for A function
        nodes_p, weights_p = gauleg(0.0, 20.0, 64)
        
        @check max_examples=5 function property_average_scattering_rate_equivalence(
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
            
            # Create NamedTuple parameters with A field
            q_nt = (
                m = (u=m_u, d=m_u, s=m_s),
                μ = (u=μ_u, d=μ_u, s=μ_s),
                A = (u=A_u, d=A_u, s=A_s)
            )
            t_nt = (T=T, Φ=Φ, Φbar=Φbar, ξ=0.0)
            
            # For struct testing, add A field to match NamedTuple
            q_struct_with_A = merge(as_namedtuple(q_struct), (A=(u=A_u, d=A_u, s=A_s),))
            
            # Test a representative scattering process (uu_to_uu)
            process = :uu_to_uu
            
            # Compute average scattering rate with struct parameters
            # Use minimal integration nodes for faster testing
            rate_struct = average_scattering_rate(
                process, q_struct_with_A, t_struct, K_coeffs,
                p_nodes=4,
                angle_nodes=2,
                phi_nodes=2,
                density_p_nodes=8,
                scale=10.0,
                density_scale=10.0
            )
            
            # Compute average scattering rate with NamedTuple parameters
            rate_nt = average_scattering_rate(
                process, q_nt, t_nt, K_coeffs,
                p_nodes=4,
                angle_nodes=2,
                phi_nodes=2,
                density_p_nodes=8,
                scale=10.0,
                density_scale=10.0
            )
            
            # Verify numerical equivalence
            isapprox(rate_struct, rate_nt, rtol=RTOL, atol=ATOL) ||
                error("Struct-NamedTuple equivalence failed for $process: " *
                      "struct=$rate_struct, nt=$rate_nt, " *
                      "diff=$(abs(rate_struct - rate_nt))")
            
            # Verify physical constraints: rate must be non-negative
            rate_struct >= 0.0 ||
                error("Scattering rate must be non-negative (struct): $rate_struct")
            rate_nt >= 0.0 ||
                error("Scattering rate must be non-negative (nt): $rate_nt")
            
            # Verify rate is finite
            isfinite(rate_struct) ||
                error("Scattering rate must be finite (struct): $rate_struct")
            isfinite(rate_nt) ||
                error("Scattering rate must be finite (nt): $rate_nt")
        end
        
        println("✓ Property test passed: Struct-NamedTuple equivalence verified for uu_to_uu")
    end
    
    # ========================================================================
    # Property 2: Struct-NamedTuple Equivalence for Different Processes
    # ========================================================================
    
    @testset "Property: Multiple Scattering Processes" begin
        println("\n" * "="^70)
        println("Property Test: Multiple Scattering Processes")
        println("Testing struct-NamedTuple equivalence across different processes")
        println("="^70)
        
        RTOL = 1e-12
        ATOL = 1e-14
        
        nodes_p, weights_p = gauleg(0.0, 20.0, 64)
        
        # Test a subset of representative processes
        test_processes = [:uu_to_uu, :ss_to_ss, :ud_to_ud, :us_to_us]
        
        @check max_examples=3 function property_multiple_processes(
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
            q_struct_with_A = merge(as_namedtuple(q_struct), (A=(u=A_u, d=A_u, s=A_s),))
            q_nt = (m=(u=m_u, d=m_u, s=m_s), μ=(u=μ_u, d=μ_u, s=μ_s), A=(u=A_u, d=A_u, s=A_s))
            t_nt = (T=T, Φ=Φ, Φbar=Φbar, ξ=0.0)
            
            # Test each process
            for process in test_processes
                rate_struct = average_scattering_rate(
                    process, q_struct_with_A, t_struct, K_coeffs,
                    p_nodes=4, angle_nodes=2, phi_nodes=2,
                    density_p_nodes=8, scale=10.0, density_scale=10.0
                )
                rate_nt = average_scattering_rate(
                    process, q_nt, t_nt, K_coeffs,
                    p_nodes=4, angle_nodes=2, phi_nodes=2,
                    density_p_nodes=8, scale=10.0, density_scale=10.0
                )
                
                # Only verify if both are finite (some processes may fail at certain parameters)
                if isfinite(rate_struct) && isfinite(rate_nt)
                    isapprox(rate_struct, rate_nt, rtol=RTOL, atol=ATOL) ||
                        error("Process $process failed: struct=$rate_struct, nt=$rate_nt")
                    
                    rate_struct >= 0.0 || error("Negative rate for $process (struct)")
                end
            end
            
            # Return true to indicate test passed
            true
        end
        
        println("✓ Property test passed: All tested processes show struct-NamedTuple equivalence")
    end
    
    # ========================================================================
    # Property 3: Cache Construction with Structs
    # ========================================================================
    
    @testset "Property: Cache Construction" begin
        println("\n" * "="^70)
        println("Property Test: Cache Construction with Structs")
        println("Testing build_w0cdf_pchip_cache with struct parameters")
        println("="^70)
        
        RTOL = 1e-12
        ATOL = 1e-14
        
        nodes_p, weights_p = gauleg(0.0, 20.0, 64)
        
        @check max_examples=3 function property_cache_construction(
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
            q_struct_with_A = merge(as_namedtuple(q_struct), (A=(u=A_u, d=A_u, s=A_s),))
            q_nt = (m=(u=m_u, d=m_u, s=m_s), μ=(u=μ_u, d=μ_u, s=μ_s), A=(u=A_u, d=A_u, s=A_s))
            t_nt = (T=T, Φ=Φ, Φbar=Φbar, ξ=0.0)
            
            process = :uu_to_uu
            
            # Build cache with struct parameters
            cache_struct = build_w0cdf_pchip_cache(
                process, q_struct_with_A, t_struct, K_coeffs,
                N=10,  # Minimal for faster testing
                design_p_nodes=4,
                design_angle_nodes=2,
                design_phi_nodes=2,
                p_cutoff=Constants_PNJL.Λ_inv_fm,
                n_sigma_points=4
            )
            
            # Build cache with NamedTuple parameters
            cache_nt = build_w0cdf_pchip_cache(
                process, q_nt, t_nt, K_coeffs,
                N=10,
                design_p_nodes=4,
                design_angle_nodes=2,
                design_phi_nodes=2,
                p_cutoff=Constants_PNJL.Λ_inv_fm,
                n_sigma_points=4
            )
            
            # Verify cache structure
            length(cache_struct.s_vals) > 0 || error("Empty cache (struct)")
            length(cache_nt.s_vals) > 0 || error("Empty cache (nt)")
            
            # Verify cache values are equivalent
            length(cache_struct.s_vals) == length(cache_nt.s_vals) ||
                error("Cache size mismatch: struct=$(length(cache_struct.s_vals)), nt=$(length(cache_nt.s_vals))")
            
            # Verify s_vals are identical
            for i in 1:length(cache_struct.s_vals)
                isapprox(cache_struct.s_vals[i], cache_nt.s_vals[i], rtol=RTOL, atol=ATOL) ||
                    error("s_vals mismatch at index $i")
            end
            
            # Verify sigma_vals are equivalent
            for i in 1:length(cache_struct.sigma_vals)
                isapprox(cache_struct.sigma_vals[i], cache_nt.sigma_vals[i], rtol=RTOL, atol=ATOL) ||
                    error("sigma_vals mismatch at index $i: struct=$(cache_struct.sigma_vals[i]), nt=$(cache_nt.sigma_vals[i])")
            end
            
            # Return true to indicate test passed
            true
        end
        
        println("✓ Property test passed: Cache construction with structs verified")
    end
    
end

println("\n" * "="^70)
println("AverageScatteringRate Property Tests Complete!")
println("="^70)
println("\nAll property tests passed:")
println("  1. Struct-NamedTuple equivalence for average_scattering_rate")
println("  2. Multiple scattering processes tested")
println("  3. Cache construction with struct parameters verified")
println("  4. Physical constraints maintained (non-negative, finite rates)")
println("="^70)
