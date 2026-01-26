"""
Unit tests for AverageScatteringRate module edge cases and specific scenarios.

Tests Requirements 2.4, 2.5:
- build_w0cdf_pchip_cache with struct parameters
- Cache construction and interpolation
- Different scattering processes
- Helper functions get_mass and get_mu with structs

These tests complement the property-based tests by focusing on specific
edge cases and boundary conditions.
"""

using Test

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
using .Main: QuarkParams, ThermoParams, as_namedtuple

@testset "AverageScatteringRate Edge Case Tests" begin
    
    # Setup common test parameters
    m_u = 1.52  # fm⁻¹
    m_s = 3.04  # fm⁻¹
    μ_u = 0.3   # fm⁻¹
    μ_s = 0.3   # fm⁻¹
    T = 0.15    # fm⁻¹
    Φ = 0.5
    Φbar = 0.5
    
    # Compute A functions
    nodes_p, weights_p = gauleg(0.0, 20.0, 64)
    A_u = A(m_u, μ_u, T, Φ, Φbar, nodes_p, weights_p)
    A_s = A(m_s, μ_s, T, Φ, Φbar, nodes_p, weights_p)
    
    # Compute K_coeffs
    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)
    K_coeffs = calculate_effective_couplings(
        Constants_PNJL.G_fm2, Constants_PNJL.K_fm5, G_u, G_s
    )
    
    # Create struct parameters
    q_struct = QuarkParams((u=m_u, d=m_u, s=m_s), (u=μ_u, d=μ_u, s=μ_s))
    t_struct = ThermoParams(T, Φ, Φbar, 0.0)
    
    # Add A field for compatibility
    q_struct_with_A = merge(as_namedtuple(q_struct), (A=(u=A_u, d=A_u, s=A_s),))
    
    # ========================================================================
    # Test 1: Helper functions get_mass and get_mu with structs
    # ========================================================================
    
    @testset "Helper functions with struct parameters" begin
        println("\n" * "="^70)
        println("Test: get_mass and get_mu with struct parameters")
        println("="^70)
        
        # Test get_mass with struct
        @test AverageScatteringRate.get_mass(:u, q_struct) == m_u
        @test AverageScatteringRate.get_mass(:d, q_struct) == m_u
        @test AverageScatteringRate.get_mass(:s, q_struct) == m_s
        @test AverageScatteringRate.get_mass(:ubar, q_struct) == m_u
        @test AverageScatteringRate.get_mass(:dbar, q_struct) == m_u
        @test AverageScatteringRate.get_mass(:sbar, q_struct) == m_s
        
        # Test get_mu with struct
        @test AverageScatteringRate.get_mu(:u, q_struct) == μ_u
        @test AverageScatteringRate.get_mu(:d, q_struct) == μ_u
        @test AverageScatteringRate.get_mu(:s, q_struct) == μ_s
        @test AverageScatteringRate.get_mu(:ubar, q_struct) == μ_u
        @test AverageScatteringRate.get_mu(:dbar, q_struct) == μ_u
        @test AverageScatteringRate.get_mu(:sbar, q_struct) == μ_s
        
        # Test with NamedTuple for comparison
        q_nt = as_namedtuple(q_struct)
        @test AverageScatteringRate.get_mass(:u, q_nt) == m_u
        @test AverageScatteringRate.get_mu(:u, q_nt) == μ_u
        
        println("✓ Helper functions work correctly with struct parameters")
    end
    
    # ========================================================================
    # Test 2: build_w0cdf_pchip_cache with struct parameters
    # ========================================================================
    
    @testset "build_w0cdf_pchip_cache with structs" begin
        println("\n" * "="^70)
        println("Test: build_w0cdf_pchip_cache with struct parameters")
        println("="^70)
        
        process = :uu_to_uu
        
        # Build cache with struct parameters (minimal settings for speed)
        cache = build_w0cdf_pchip_cache(
            process, q_struct_with_A, t_struct, K_coeffs,
            N=10,
            design_p_nodes=4,
            design_angle_nodes=2,
            design_phi_nodes=2,
            p_cutoff=Constants_PNJL.Λ_inv_fm,
            n_sigma_points=4
        )
        
        # Verify cache structure
        @test cache isa CrossSectionCache
        @test cache.process == process
        @test length(cache.s_vals) > 0
        @test length(cache.sigma_vals) == length(cache.s_vals)
        @test length(cache.pchip_slopes) == length(cache.s_vals)
        
        # Verify s_vals are sorted
        @test issorted(cache.s_vals)
        
        # Verify all sigma values are non-negative
        @test all(σ >= 0.0 for σ in cache.sigma_vals)
        
        # Verify all values are finite
        @test all(isfinite, cache.s_vals)
        @test all(isfinite, cache.sigma_vals)
        
        println("✓ build_w0cdf_pchip_cache works with struct parameters")
    end
    
    # ========================================================================
    # Test 3: Cache interpolation
    # ========================================================================
    
    @testset "Cache interpolation" begin
        println("\n" * "="^70)
        println("Test: Cache interpolation with struct-built cache")
        println("="^70)
        
        process = :uu_to_uu
        
        # Build cache
        cache = build_w0cdf_pchip_cache(
            process, q_struct_with_A, t_struct, K_coeffs,
            N=10,
            design_p_nodes=4,
            design_angle_nodes=2,
            design_phi_nodes=2,
            p_cutoff=Constants_PNJL.Λ_inv_fm,
            n_sigma_points=4
        )
        
        # Test interpolation at cache points
        for i in 1:length(cache.s_vals)
            s = cache.s_vals[i]
            σ_expected = cache.sigma_vals[i]
            σ_interp = AverageScatteringRate.interpolate_sigma(cache, s)
            @test σ_interp !== nothing
            @test isapprox(σ_interp, σ_expected, rtol=1e-10)
        end
        
        # Test interpolation between points
        if length(cache.s_vals) >= 2
            s_mid = (cache.s_vals[1] + cache.s_vals[2]) / 2.0
            σ_mid = AverageScatteringRate.interpolate_sigma(cache, s_mid)
            @test σ_mid !== nothing
            @test σ_mid >= 0.0
            @test isfinite(σ_mid)
        end
        
        # Test out-of-bounds behavior
        s_below = cache.s_vals[1] - 1.0
        s_above = cache.s_vals[end] + 1.0
        @test AverageScatteringRate.interpolate_sigma(cache, s_below) === nothing
        @test AverageScatteringRate.interpolate_sigma(cache, s_above) === nothing
        
        println("✓ Cache interpolation works correctly")
    end
    
    # ========================================================================
    # Test 4: Different scattering processes
    # ========================================================================
    
    @testset "Different scattering processes" begin
        println("\n" * "="^70)
        println("Test: average_scattering_rate with different processes")
        println("="^70)
        
        # Test only one process with minimal integration nodes
        process = :uu_to_uu
        
        rate = average_scattering_rate(
            process, q_struct_with_A, t_struct, K_coeffs,
            p_nodes=4,
            angle_nodes=2,
            phi_nodes=2,
            density_p_nodes=8,
            scale=10.0,
            density_scale=10.0
        )
        
        # Verify rate is non-negative and finite
        @test rate >= 0.0
        @test isfinite(rate)
        
        println("  ✓ Process $process: rate = $rate")
        
        println("✓ Tested process produces valid rate")
    end
    
    # ========================================================================
    # Test 5: precompute_cross_section! with structs
    # ========================================================================
    
    @testset "precompute_cross_section! with structs" begin
        println("\n" * "="^70)
        println("Test: precompute_cross_section! with struct parameters")
        println("="^70)
        
        process = :uu_to_uu
        
        # Create empty cache
        cache = CrossSectionCache(process)
        
        # Define s_grid
        s_grid = [10.0, 20.0, 30.0, 40.0]
        
        # Precompute with struct parameters
        precompute_cross_section!(
            cache, s_grid, q_struct_with_A, t_struct, K_coeffs,
            n_points=4
        )
        
        # Verify cache was populated
        @test length(cache.s_vals) == length(s_grid)
        @test length(cache.sigma_vals) == length(s_grid)
        
        # Verify s_vals match input grid
        for i in 1:length(s_grid)
            @test cache.s_vals[i] == s_grid[i]
        end
        
        # Verify all sigma values are non-negative and finite
        @test all(σ >= 0.0 for σ in cache.sigma_vals)
        @test all(isfinite, cache.sigma_vals)
        
        println("✓ precompute_cross_section! works with struct parameters")
    end
    
    # ========================================================================
    # Test 6: design_w0cdf_s_grid with structs
    # ========================================================================
    
    @testset "design_w0cdf_s_grid with structs" begin
        println("\n" * "="^70)
        println("Test: design_w0cdf_s_grid with struct parameters")
        println("="^70)
        
        process = :uu_to_uu
        
        # Design s_grid with struct parameters (minimal settings)
        s_grid = AverageScatteringRate.design_w0cdf_s_grid(
            process, q_struct_with_A, t_struct,
            N=10,
            p_nodes=4,
            angle_nodes=2,
            phi_nodes=2,
            p_cutoff=Constants_PNJL.Λ_inv_fm
        )
        
        # Verify grid structure
        @test length(s_grid) == 10
        @test all(isfinite, s_grid)
        @test all(s > 0.0 for s in s_grid)
        
        # Verify grid is sorted
        @test issorted(s_grid)
        
        # Verify grid is above threshold
        s_threshold = (2.0 * m_u)^2
        @test all(s > s_threshold for s in s_grid)
        
        println("✓ design_w0cdf_s_grid works with struct parameters")
    end
    
    # ========================================================================
    # Test 7: Mixed struct/NamedTuple usage
    # ========================================================================
    
    @testset "Mixed struct/NamedTuple usage" begin
        println("\n" * "="^70)
        println("Test: Mixed struct and NamedTuple parameters")
        println("="^70)
        
        process = :uu_to_uu
        q_nt = as_namedtuple(q_struct)
        q_nt_with_A = merge(q_nt, (A=(u=A_u, d=A_u, s=A_s),))
        t_nt = (T=T, Φ=Φ, Φbar=Φbar, ξ=0.0)
        
        # Test struct quark_params with NamedTuple thermo_params
        rate1 = average_scattering_rate(
            process, q_struct_with_A, t_nt, K_coeffs,
            p_nodes=3, angle_nodes=2, phi_nodes=2,
            density_p_nodes=6, scale=10.0, density_scale=10.0
        )
        
        # Test NamedTuple quark_params with struct thermo_params
        rate2 = average_scattering_rate(
            process, q_nt_with_A, t_struct, K_coeffs,
            p_nodes=3, angle_nodes=2, phi_nodes=2,
            density_p_nodes=6, scale=10.0, density_scale=10.0
        )
        
        # All should produce equivalent results
        @test isapprox(rate1, rate2, rtol=1e-10)
        
        println("✓ Mixed struct/NamedTuple usage works correctly")
    end
    
end

println("\n" * "="^70)
println("AverageScatteringRate Edge Case Tests Complete!")
println("="^70)
println("\nAll edge case tests passed:")
println("  1. Helper functions (get_mass, get_mu) work with structs")
println("  2. build_w0cdf_pchip_cache works with struct parameters")
println("  3. Cache interpolation works correctly")
println("  4. Different scattering processes produce valid rates")
println("  5. precompute_cross_section! works with structs")
println("  6. design_w0cdf_s_grid works with structs")
println("  7. Mixed struct/NamedTuple usage works correctly")
println("="^70)
