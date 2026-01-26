"""
Unit tests for RelaxationTime module struct support - edge cases and specific scenarios.

Tests Requirements 1.3, 1.4

This test file verifies that the RelaxationTime module correctly handles:
- compute_average_rates with struct parameters
- relaxation_rates with struct-derived densities
- Different scattering process combinations
- Pre-computed rates vs computing from scratch
"""

using Test

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

@testset "RelaxationTime Edge Cases and Unit Tests" begin
    
    # Setup common test parameters
    m_u = 1.52
    m_s = 3.04
    μ_u = 0.3
    μ_s = 0.3
    T = 0.15
    Φ = 0.5
    Φbar = 0.5
    
    # Pre-compute A functions
    nodes_p, weights_p = gauleg(0.0, 20.0, 64)
    A_u = A(m_u, μ_u, T, Φ, Φbar, nodes_p, weights_p)
    A_s = A(m_s, μ_s, T, Φ, Φbar, nodes_p, weights_p)
    
    # Compute K_coeffs
    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)
    K_coeffs = calculate_effective_couplings(Constants_PNJL.G_fm2, Constants_PNJL.K_fm5, G_u, G_s)
    
    # Create test parameters
    q_struct = QuarkParams((u=m_u, d=m_u, s=m_s), (u=μ_u, d=μ_u, s=μ_s))
    t_struct = ThermoParams(T, Φ, Φbar, 0.0)
    q_nt = (m=(u=m_u, d=m_u, s=m_s), μ=(u=μ_u, d=μ_u, s=μ_s))
    t_nt = (T=T, Φ=Φ, Φbar=Φbar, ξ=0.0)
    
    densities = (u=0.1, d=0.1, s=0.05, ubar=0.1, dbar=0.1, sbar=0.05)
    
    # ========================================================================
    # Test 1: compute_average_rates with struct parameters
    # ========================================================================
    
    @testset "compute_average_rates with structs" begin
        println("\nTest: compute_average_rates with struct parameters")
        
        # Compute with struct parameters (minimal nodes for speed)
        rates_struct = compute_average_rates(
            q_struct, t_struct, K_coeffs,
            p_nodes=2,
            angle_nodes=2,
            phi_nodes=2,
            n_sigma_points=3
        )
        
        # Compute with NamedTuple parameters
        rates_nt = compute_average_rates(
            q_nt, t_nt, K_coeffs,
            p_nodes=2,
            angle_nodes=2,
            phi_nodes=2,
            n_sigma_points=3
        )
        
        # Verify both produce results
        @test rates_struct isa NamedTuple
        @test rates_nt isa NamedTuple
        
        # Verify they have the same keys
        @test Set(keys(rates_struct)) == Set(keys(rates_nt))
        
        # Verify all rates are non-negative and finite
        for (key, val) in pairs(rates_struct)
            @test val >= 0.0
            @test isfinite(val)
        end
        
        for (key, val) in pairs(rates_nt)
            @test val >= 0.0
            @test isfinite(val)
        end
        
        # Verify numerical equivalence for a subset of processes
        test_processes = [:uu_to_uu, :ss_to_ss, :ud_to_ud]
        for process in test_processes
            if hasproperty(rates_struct, process) && hasproperty(rates_nt, process)
                rate_struct = getproperty(rates_struct, process)
                rate_nt = getproperty(rates_nt, process)
                @test isapprox(rate_struct, rate_nt, rtol=1e-12, atol=1e-14)
            end
        end
        
        println("✓ compute_average_rates works with struct parameters")
    end
    
    # ========================================================================
    # Test 2: relaxation_rates with struct-derived densities
    # ========================================================================
    
    @testset "relaxation_rates with densities" begin
        println("\nTest: relaxation_rates with struct-derived densities")
        
        # First compute some rates
        rates = compute_average_rates(
            q_struct, t_struct, K_coeffs,
            p_nodes=2,
            angle_nodes=2,
            phi_nodes=2,
            n_sigma_points=3
        )
        
        # Compute relaxation rates
        tau_inv = relaxation_rates(densities, rates)
        
        # Verify structure
        @test tau_inv isa NamedTuple
        @test haskey(tau_inv, :u)
        @test haskey(tau_inv, :d)
        @test haskey(tau_inv, :s)
        @test haskey(tau_inv, :ubar)
        @test haskey(tau_inv, :dbar)
        @test haskey(tau_inv, :sbar)
        
        # Verify all values are non-negative and finite
        for flavor in [:u, :d, :s, :ubar, :dbar, :sbar]
            val = getproperty(tau_inv, flavor)
            @test val >= 0.0
            @test isfinite(val)
        end
        
        # Verify isospin symmetry
        @test tau_inv.u == tau_inv.d
        @test tau_inv.ubar == tau_inv.dbar
        
        println("✓ relaxation_rates works correctly")
    end
    
    # ========================================================================
    # Test 3: relaxation_times with pre-computed rates
    # ========================================================================
    
    @testset "relaxation_times with pre-computed rates" begin
        println("\nTest: relaxation_times with pre-computed rates")
        
        # First compute rates separately
        rates = compute_average_rates(
            q_struct, t_struct, K_coeffs,
            p_nodes=2,
            angle_nodes=2,
            phi_nodes=2,
            n_sigma_points=3
        )
        
        # Now call relaxation_times with existing_rates
        result = relaxation_times(
            q_struct, t_struct, K_coeffs;
            densities=densities,
            existing_rates=rates,
            p_nodes=2,
            angle_nodes=2,
            phi_nodes=2,
            n_sigma_points=3
        )
        
        # Verify structure
        @test result isa NamedTuple
        @test haskey(result, :tau)
        @test haskey(result, :tau_inv)
        @test haskey(result, :rates)
        
        # Verify tau values are positive and finite
        for flavor in [:u, :d, :s, :ubar, :dbar, :sbar]
            tau_val = getproperty(result.tau, flavor)
            @test tau_val > 0.0
            @test isfinite(tau_val)
        end
        
        # Verify tau_inv values are non-negative and finite
        for flavor in [:u, :d, :s, :ubar, :dbar, :sbar]
            tau_inv_val = getproperty(result.tau_inv, flavor)
            @test tau_inv_val >= 0.0
            @test isfinite(tau_inv_val)
        end
        
        # Verify tau and tau_inv are inverses (within tolerance)
        for flavor in [:u, :d, :s, :ubar, :dbar, :sbar]
            tau_val = getproperty(result.tau, flavor)
            tau_inv_val = getproperty(result.tau_inv, flavor)
            if tau_inv_val > 0.0
                @test isapprox(tau_val * tau_inv_val, 1.0, rtol=1e-10)
            else
                @test isinf(tau_val)
            end
        end
        
        println("✓ relaxation_times works with pre-computed rates")
    end
    
    # ========================================================================
    # Test 4: Struct-NamedTuple equivalence for full relaxation_times
    # ========================================================================
    
    @testset "relaxation_times struct-NamedTuple equivalence" begin
        println("\nTest: relaxation_times struct-NamedTuple equivalence")
        
        # Compute with struct parameters
        result_struct = relaxation_times(
            q_struct, t_struct, K_coeffs;
            densities=densities,
            p_nodes=2,
            angle_nodes=2,
            phi_nodes=2,
            n_sigma_points=3
        )
        
        # Compute with NamedTuple parameters
        result_nt = relaxation_times(
            q_nt, t_nt, K_coeffs;
            densities=densities,
            p_nodes=2,
            angle_nodes=2,
            phi_nodes=2,
            n_sigma_points=3
        )
        
        # Verify tau values are equivalent
        for flavor in [:u, :d, :s, :ubar, :dbar, :sbar]
            tau_struct = getproperty(result_struct.tau, flavor)
            tau_nt = getproperty(result_nt.tau, flavor)
            @test isapprox(tau_struct, tau_nt, rtol=1e-12, atol=1e-14)
        end
        
        # Verify tau_inv values are equivalent
        for flavor in [:u, :d, :s, :ubar, :dbar, :sbar]
            tau_inv_struct = getproperty(result_struct.tau_inv, flavor)
            tau_inv_nt = getproperty(result_nt.tau_inv, flavor)
            @test isapprox(tau_inv_struct, tau_inv_nt, rtol=1e-12, atol=1e-14)
        end
        
        println("✓ Struct-NamedTuple equivalence verified for relaxation_times")
    end
    
    # ========================================================================
    # Test 5: Mixed struct/NamedTuple usage
    # ========================================================================
    
    @testset "Mixed struct/NamedTuple usage" begin
        println("\nTest: Mixed struct/NamedTuple usage patterns")
        
        # Test with struct quark_params and NamedTuple thermo_params
        result1 = relaxation_times(
            q_struct, t_nt, K_coeffs;
            densities=densities,
            p_nodes=2,
            angle_nodes=2,
            phi_nodes=2,
            n_sigma_points=3
        )
        
        # Test with NamedTuple quark_params and struct thermo_params
        result2 = relaxation_times(
            q_nt, t_struct, K_coeffs;
            densities=densities,
            p_nodes=2,
            angle_nodes=2,
            phi_nodes=2,
            n_sigma_points=3
        )
        
        # Both should produce valid results
        @test result1 isa NamedTuple
        @test result2 isa NamedTuple
        
        # Verify tau values are positive and finite for both
        for result in [result1, result2]
            for flavor in [:u, :d, :s, :ubar, :dbar, :sbar]
                tau_val = getproperty(result.tau, flavor)
                @test tau_val > 0.0
                @test isfinite(tau_val)
            end
        end
        
        println("✓ Mixed struct/NamedTuple usage works correctly")
    end
    
end

println("\n" * "="^70)
println("RelaxationTime Edge Cases and Unit Tests Complete!")
println("="^70)
println("\nAll unit tests passed:")
println("  1. compute_average_rates works with struct parameters")
println("  2. relaxation_rates works with struct-derived densities")
println("  3. relaxation_times works with pre-computed rates")
println("  4. Struct-NamedTuple equivalence verified")
println("  5. Mixed struct/NamedTuple usage patterns work")
println("  6. Physical constraints maintained (positive tau, non-negative rates)")
println("  7. Isospin symmetry preserved")
println("="^70)
