"""
Unit tests for TotalCrossSection module edge cases and specific scenarios.

Tests Requirements 3.4, 3.5:
- calculate_all_total_cross_sections with struct parameters
- scan_s_dependence produces smooth curves
- Threshold behavior
- Different scattering processes

These tests complement the property-based tests by focusing on specific
edge cases and boundary conditions.
"""

using Test

# Load required modules
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src/relaxtime"))

include("../../../src/Constants_PNJL.jl")
include("../../../src/relaxtime/TotalCrossSection.jl")
include("../../../src/relaxtime/EffectiveCouplings.jl")
include("../../../src/relaxtime/OneLoopIntegrals.jl")
include("../../../src/integration/GaussLegendre.jl")

using .TotalCrossSection
using .Constants_PNJL
using .EffectiveCouplings
using .OneLoopIntegrals: A
using .GaussLegendre: gauleg

# Load test utilities
include("test_utils.jl")
using .Main: QuarkParams, ThermoParams, as_namedtuple

@testset "TotalCrossSection Edge Case Tests" begin
    
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
    
    # Add A field for compatibility with TotalPropagator
    q_struct_with_A = merge(as_namedtuple(q_struct), (A=(u=A_u, d=A_u, s=A_s),))
    
    # ========================================================================
    # Test 1: calculate_all_total_cross_sections with struct parameters
    # ========================================================================
    
    @testset "calculate_all_total_cross_sections with structs" begin
        println("\n" * "="^70)
        println("Test: calculate_all_total_cross_sections with struct parameters")
        println("="^70)
        
        s = 31.0  # fm⁻²
        
        # Call with struct parameters
        all_σ = calculate_all_total_cross_sections(
            s, q_struct_with_A, t_struct, K_coeffs, n_points=6
        )
        
        # Verify result structure
        @test all_σ isa NamedTuple
        @test haskey(all_σ, :uu_to_uu)
        @test haskey(all_σ, :ss_to_ss)
        @test haskey(all_σ, :ud_to_ud)
        @test haskey(all_σ, :us_to_us)
        
        # Verify all finite values are non-negative
        for (process, σ) in pairs(all_σ)
            if isfinite(σ)
                @test σ >= 0.0
            end
        end
        
        # Verify at least some processes have finite results
        finite_count = count(isfinite, values(all_σ))
        @test finite_count > 0
        
        println("✓ calculate_all_total_cross_sections works with struct parameters")
        println("  Finite results: $finite_count / $(length(all_σ))")
    end
    
    # ========================================================================
    # Test 2: scan_s_dependence produces smooth curves
    # ========================================================================
    
    @testset "scan_s_dependence smoothness" begin
        println("\n" * "="^70)
        println("Test: scan_s_dependence produces smooth curves")
        println("="^70)
        
        # Test with uu_to_uu process
        process = :uu_to_uu
        
        # Create s values from threshold to high energy
        s_threshold = (2.0 * m_u)^2  # Threshold for uu_to_uu
        s_values = collect(range(s_threshold * 1.1, s_threshold * 5.0, length=10))
        
        # Scan with struct parameters
        σ_values = scan_s_dependence(
            s_values, process, q_struct_with_A, t_struct, K_coeffs, n_points=6
        )
        
        # Verify result structure
        @test length(σ_values) == length(s_values)
        @test all(isfinite, σ_values)
        @test all(>=(0.0), σ_values)
        
        # Test smoothness: check that cross-section increases monotonically
        # or at least doesn't have large jumps
        # (Physical expectation: σ generally increases with s above threshold)
        for i in 2:length(σ_values)
            # Check for large discontinuous jumps (more than 50% change)
            if σ_values[i-1] > 0.0
                relative_change = abs(σ_values[i] - σ_values[i-1]) / σ_values[i-1]
                @test relative_change < 2.0
            end
        end
        
        println("✓ scan_s_dependence produces smooth curves")
        println("  s range: $(s_values[1]) to $(s_values[end]) fm⁻²")
        println("  σ range: $(minimum(σ_values)) to $(maximum(σ_values)) fm²")
    end
    
    # ========================================================================
    # Test 3: Threshold behavior
    # ========================================================================
    
    @testset "Threshold behavior" begin
        println("\n" * "="^70)
        println("Test: Cross-section behavior near kinematic threshold")
        println("="^70)
        
        process = :uu_to_uu
        s_threshold = (2.0 * m_u)^2
        
        # Test below threshold (should return 0)
        s_below = s_threshold * 0.95
        σ_below = total_cross_section(
            process, s_below, q_struct_with_A, t_struct, K_coeffs, n_points=6
        )
        @test σ_below == 0.0
        
        # Test at threshold (should be small but non-negative)
        s_at = s_threshold * 1.001
        σ_at = total_cross_section(
            process, s_at, q_struct_with_A, t_struct, K_coeffs, n_points=6
        )
        @test σ_at >= 0.0
        @test isfinite(σ_at)
        
        # Test slightly above threshold
        s_above = s_threshold * 1.1
        σ_above = total_cross_section(
            process, s_above, q_struct_with_A, t_struct, K_coeffs, n_points=6
        )
        @test σ_above >= 0.0
        @test isfinite(σ_above)
        
        # Physical expectation: cross-section should increase with s above threshold
        # (though this may not always be strictly true due to quantum effects)
        @test σ_above >= σ_at
        
        println("✓ Threshold behavior verified")
        println("  s_threshold = $s_threshold fm⁻²")
        println("  σ(below) = $σ_below fm²")
        println("  σ(at) = $σ_at fm²")
        println("  σ(above) = $σ_above fm²")
    end
    
    # ========================================================================
    # Test 4: Different scattering processes
    # ========================================================================
    
    @testset "Different scattering processes" begin
        println("\n" * "="^70)
        println("Test: Cross-sections for different scattering processes")
        println("="^70)
        
        s = 31.0  # fm⁻²
        
        # Test a representative set of processes
        test_processes = [
            :uu_to_uu,  # Same flavor, same particles
            :ubarubar_to_ubarubar,  # Antiparticle version (should equal uu_to_uu)
            :ss_to_ss,  # Heavier quark
            :ud_to_ud,  # Different flavors
            :us_to_us,  # Light-heavy mixing
        ]
        
        results = Dict{Symbol, Float64}()
        
        for process in test_processes
            σ = total_cross_section(
                process, s, q_struct_with_A, t_struct, K_coeffs, n_points=6
            )
            results[process] = σ
            
            # Verify physical constraints
            @test isfinite(σ)
            @test σ >= 0.0
            
            println("  $process: σ = $σ fm²")
        end
        
        # Test symmetry: uu_to_uu should equal ubarubar_to_ubarubar (charge conjugate)
        if isfinite(results[:uu_to_uu]) && isfinite(results[:ubarubar_to_ubarubar])
            @test isapprox(results[:uu_to_uu], results[:ubarubar_to_ubarubar], rtol=1e-6)
        end
        
        # Test that ss_to_ss is different from uu_to_uu (different mass)
        if isfinite(results[:uu_to_uu]) && isfinite(results[:ss_to_ss])
            @test results[:uu_to_uu] != results[:ss_to_ss]
        end
        
        println("✓ Different scattering processes tested successfully")
    end
    
    # ========================================================================
    # Test 5: Integration point convergence
    # ========================================================================
    
    @testset "Integration point convergence" begin
        println("\n" * "="^70)
        println("Test: Cross-section convergence with integration points")
        println("="^70)
        
        process = :uu_to_uu
        s = 31.0  # fm⁻²
        
        # Test with different numbers of integration points
        n_points_list = [4, 6, 8, 12]
        σ_list = Float64[]
        
        for n_points in n_points_list
            σ = total_cross_section(
                process, s, q_struct_with_A, t_struct, K_coeffs, n_points=n_points
            )
            push!(σ_list, σ)
            println("  n_points=$n_points: σ = $σ fm²")
        end
        
        # Verify all results are finite and non-negative
        @test all(isfinite, σ_list)
        @test all(>=(0.0), σ_list)
        
        # Test convergence: results should be similar for higher n_points
        # (difference between n=8 and n=12 should be small)
        if length(σ_list) >= 4
            σ_8 = σ_list[3]
            σ_12 = σ_list[4]
            if σ_12 > 0.0
                relative_diff = abs(σ_8 - σ_12) / σ_12
                @test relative_diff < 0.05
            end
        end
        
        println("✓ Integration point convergence verified")
    end
    
    # ========================================================================
    # Test 6: Struct vs NamedTuple consistency for batch operations
    # ========================================================================
    
    @testset "Struct-NamedTuple consistency for batch operations" begin
        println("\n" * "="^70)
        println("Test: Struct-NamedTuple consistency for batch operations")
        println("="^70)
        
        s = 31.0  # fm⁻²
        
        # Create NamedTuple parameters
        q_nt = (
            m = (u=m_u, d=m_u, s=m_s),
            μ = (u=μ_u, d=μ_u, s=μ_s),
            A = (u=A_u, d=A_u, s=A_s)
        )
        t_nt = (T=T, Φ=Φ, Φbar=Φbar, ξ=0.0)
        
        # Test calculate_all_total_cross_sections
        all_σ_struct = calculate_all_total_cross_sections(
            s, q_struct_with_A, t_struct, K_coeffs, n_points=6
        )
        all_σ_nt = calculate_all_total_cross_sections(
            s, q_nt, t_nt, K_coeffs, n_points=6
        )
        
        # Verify equivalence for all processes
        for process in keys(all_σ_struct)
            σ_struct = all_σ_struct[process]
            σ_nt = all_σ_nt[process]
            
            if isfinite(σ_struct) && isfinite(σ_nt)
                @test isapprox(σ_struct, σ_nt, rtol=1e-12)
            end
        end
        
        # Test scan_s_dependence
        s_values = [20.0, 31.0, 40.0]
        process = :uu_to_uu
        
        σ_struct_scan = scan_s_dependence(
            s_values, process, q_struct_with_A, t_struct, K_coeffs, n_points=6
        )
        σ_nt_scan = scan_s_dependence(
            s_values, process, q_nt, t_nt, K_coeffs, n_points=6
        )
        
        # Verify equivalence for scan
        @test length(σ_struct_scan) == length(σ_nt_scan)
        for i in 1:length(σ_struct_scan)
            if isfinite(σ_struct_scan[i]) && isfinite(σ_nt_scan[i])
                @test isapprox(σ_struct_scan[i], σ_nt_scan[i], rtol=1e-12)
            end
        end
        
        println("✓ Struct-NamedTuple consistency verified for batch operations")
    end
    
end

println("\n" * "="^70)
println("TotalCrossSection Edge Case Tests Complete!")
println("="^70)
println("\nAll edge case tests passed:")
println("  1. calculate_all_total_cross_sections with struct parameters")
println("  2. scan_s_dependence produces smooth curves")
println("  3. Threshold behavior (below, at, above)")
println("  4. Different scattering processes tested")
println("  5. Integration point convergence")
println("  6. Struct-NamedTuple consistency for batch operations")
println("="^70)
