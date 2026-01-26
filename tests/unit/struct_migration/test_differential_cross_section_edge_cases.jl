"""
Unit tests for DifferentialCrossSection module edge cases.

Tests edge cases for:
- Threshold behavior (s near kinematic threshold)
- High-energy behavior
- Different scattering processes
- Requirements: 5.1, 5.5
"""

using Test

# Load required modules
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src/relaxtime"))

include("../../../src/relaxtime/DifferentialCrossSection.jl")

using .DifferentialCrossSection

# Load test utilities
include("test_utils.jl")
using .Main: QuarkParams, ThermoParams, as_namedtuple, approx_equal

@testset "DifferentialCrossSection Edge Cases" begin
    
    # ========================================================================
    # Test 1: Threshold Behavior (s near kinematic threshold)
    # ========================================================================
    @testset "Threshold behavior" begin
        println("\n" * "="^70)
        println("Testing threshold behavior (s near kinematic threshold)")
        println("="^70)
        
        # Test parameters
        m1 = 1.52  # u quark mass [fm⁻¹]
        m2 = 1.52  # u quark mass [fm⁻¹]
        s_threshold = (m1 + m2)^2  # Threshold: s = (m1 + m2)²
        
        # Test exactly at threshold
        @testset "Exactly at threshold" begin
            s = s_threshold
            result = check_kinematic_threshold(s, m1, m2, warn_close=false)
            @test result == true
            
            # Calculate s_12_plus and s_12_minus at threshold
            s_12_plus = s - (m1 + m2)^2  # Should be 0
            s_12_minus = s - (m1 - m2)^2  # Should be s (since m1 = m2)
            
            @test s_12_plus ≈ 0.0 atol=1e-14
            @test s_12_minus ≈ s
        end
        
        # Test just above threshold
        @testset "Just above threshold" begin
            # Test at various small distances above threshold
            for delta in [1e-6, 1e-4, 1e-2, 0.1, 1.0]
                s = s_threshold + delta
                result = check_kinematic_threshold(s, m1, m2, warn_close=false)
                @test result == true
                
                s_12_plus = s - (m1 + m2)^2
                s_12_minus = s - (m1 - m2)^2
                
                @test s_12_plus ≈ delta rtol=1e-6
                @test s_12_plus > 0.0
                @test s_12_minus > 0.0
                
                # Test that differential cross-section is finite and positive
                M_squared = 1.0  # Arbitrary matrix element
                dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
                @test isfinite(dsigma_dt)
                @test dsigma_dt > 0.0
                
                # Verify formula: dσ/dt = M²/(16π s_12_plus s_12_minus)
                expected = M_squared / (16π * s_12_plus * s_12_minus)
                @test dsigma_dt ≈ expected rtol=1e-12
            end
        end
        
        # Test below threshold (should fail)
        @testset "Below threshold" begin
            s = s_threshold - 0.1
            result = check_kinematic_threshold(s, m1, m2, warn_close=false)
            @test result == false
            
            # Attempting to calculate cross-section below threshold should error
            s_12_plus = s - (m1 + m2)^2  # Will be negative
            s_12_minus = s - (m1 - m2)^2
            M_squared = 1.0
            
            @test s_12_plus < 0.0
            @test_throws ErrorException differential_cross_section(s_12_plus, s_12_minus, M_squared)
        end
        
        # Test with different mass combinations
        @testset "Threshold with different masses" begin
            # u-s scattering (different masses)
            m1_us = 1.52  # u quark
            m2_us = 3.04  # s quark
            s_threshold_us = (m1_us + m2_us)^2
            
            # At threshold
            s = s_threshold_us
            result = check_kinematic_threshold(s, m1_us, m2_us, warn_close=false)
            @test result == true
            
            s_12_plus = s - (m1_us + m2_us)^2
            s_12_minus = s - (m1_us - m2_us)^2
            
            @test s_12_plus ≈ 0.0 atol=1e-14
            @test s_12_minus > 0.0  # Non-zero since m1 ≠ m2
            
            # Just above threshold
            s = s_threshold_us + 0.01
            s_12_plus = s - (m1_us + m2_us)^2
            s_12_minus = s - (m1_us - m2_us)^2
            M_squared = 1.0
            
            dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
            @test isfinite(dsigma_dt)
            @test dsigma_dt > 0.0
        end
        
        # Test near-threshold divergence behavior
        @testset "Near-threshold divergence" begin
            # As s → s_threshold, s_12_plus → 0, so dσ/dt → ∞
            # Test that cross-section increases as we approach threshold
            
            m1 = 1.52
            m2 = 1.52
            s_threshold = (m1 + m2)^2
            M_squared = 1.0
            
            # Calculate cross-sections at different distances from threshold
            deltas = [1.0, 0.1, 0.01, 0.001]
            cross_sections = Float64[]
            
            for delta in deltas
                s = s_threshold + delta
                s_12_plus = s - (m1 + m2)^2
                s_12_minus = s - (m1 - m2)^2
                dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
                push!(cross_sections, dsigma_dt)
            end
            
            # Verify that cross-section increases as we approach threshold
            for i in 1:(length(cross_sections)-1)
                @test cross_sections[i] < cross_sections[i+1]
            end
            
            # Verify scaling: dσ/dt ∝ 1/s_12_plus, so if s_12_plus changes by factor k,
            # cross-section changes by factor 1/k
            # cross_sections[2] / cross_sections[1] ≈ (s_12_plus_1 / s_12_plus_2)
            # Since s_12_plus_i = delta_i, we have:
            # cross_sections[2] / cross_sections[1] ≈ deltas[1] / deltas[2]
            @test cross_sections[2] / cross_sections[1] ≈ deltas[1] / deltas[2] rtol=0.1
        end
        
        println("✓ Threshold behavior tests passed")
    end
    
    # ========================================================================
    # Test 2: High-Energy Behavior
    # ========================================================================
    @testset "High-energy behavior" begin
        println("\n" * "="^70)
        println("Testing high-energy behavior")
        println("="^70)
        
        # Test parameters
        m1 = 1.52  # u quark mass [fm⁻¹]
        m2 = 1.52  # u quark mass [fm⁻¹]
        
        @testset "High-energy regime" begin
            # Test at various high energies
            for s in [100.0, 500.0, 1000.0, 5000.0, 10000.0]
                s_12_plus = s - (m1 + m2)^2
                s_12_minus = s - (m1 - m2)^2
                M_squared = 1.0
                
                # Check threshold condition
                @test check_kinematic_threshold(s, m1, m2, warn_close=false) == true
                
                # Calculate cross-section
                dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
                
                # Verify physical properties
                @test isfinite(dsigma_dt)
                @test dsigma_dt > 0.0
                
                # Verify formula
                expected = M_squared / (16π * s_12_plus * s_12_minus)
                @test dsigma_dt ≈ expected rtol=1e-12
                
                # At high energy, s_12_plus ≈ s_12_minus ≈ s (masses negligible)
                @test s_12_plus ≈ s atol=10.0  # Within 10 fm⁻² for s >> m²
                @test s_12_minus ≈ s atol=10.0
            end
        end
        
        @testset "High-energy scaling" begin
            # Test that cross-section scales correctly with energy
            m1 = 1.52
            m2 = 1.52
            M_squared = 1.0
            
            s1 = 1000.0
            s2 = 2000.0
            
            s_12_plus_1 = s1 - (m1 + m2)^2
            s_12_minus_1 = s1 - (m1 - m2)^2
            dsigma_1 = differential_cross_section(s_12_plus_1, s_12_minus_1, M_squared)
            
            s_12_plus_2 = s2 - (m1 + m2)^2
            s_12_minus_2 = s2 - (m1 - m2)^2
            dsigma_2 = differential_cross_section(s_12_plus_2, s_12_minus_2, M_squared)
            
            # At high energy, dσ/dt ∝ 1/s²
            # So dsigma_1 / dsigma_2 ≈ (s2/s1)²
            ratio = dsigma_1 / dsigma_2
            expected_ratio = (s_12_plus_2 * s_12_minus_2) / (s_12_plus_1 * s_12_minus_1)
            @test ratio ≈ expected_ratio rtol=1e-10
        end
        
        @testset "High-energy with different masses" begin
            # Test u-s scattering at high energy
            m1 = 1.52  # u quark
            m2 = 3.04  # s quark
            
            for s in [100.0, 1000.0, 10000.0]
                s_12_plus = s - (m1 + m2)^2
                s_12_minus = s - (m1 - m2)^2
                M_squared = 1.0
                
                dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
                
                @test isfinite(dsigma_dt)
                @test dsigma_dt > 0.0
                
                # Verify formula
                expected = M_squared / (16π * s_12_plus * s_12_minus)
                @test dsigma_dt ≈ expected rtol=1e-12
            end
        end
        
        @testset "Matrix element scaling at high energy" begin
            # Test that cross-section scales linearly with M²
            m1 = 1.52
            m2 = 1.52
            s = 1000.0
            
            s_12_plus = s - (m1 + m2)^2
            s_12_minus = s - (m1 - m2)^2
            
            M_squared_1 = 1.0
            M_squared_2 = 10.0
            M_squared_3 = 100.0
            
            dsigma_1 = differential_cross_section(s_12_plus, s_12_minus, M_squared_1)
            dsigma_2 = differential_cross_section(s_12_plus, s_12_minus, M_squared_2)
            dsigma_3 = differential_cross_section(s_12_plus, s_12_minus, M_squared_3)
            
            # Verify linear scaling with M²
            @test dsigma_2 / dsigma_1 ≈ M_squared_2 / M_squared_1 rtol=1e-12
            @test dsigma_3 / dsigma_1 ≈ M_squared_3 / M_squared_1 rtol=1e-12
        end
        
        println("✓ High-energy behavior tests passed")
    end
    
    # ========================================================================
    # Test 3: Different Scattering Processes
    # ========================================================================
    @testset "Different scattering processes" begin
        println("\n" * "="^70)
        println("Testing different scattering processes")
        println("="^70)
        
        # Create test parameters
        q_struct = QuarkParams((u=1.52, d=1.52, s=3.04), (u=0.3, d=0.3, s=0.3))
        
        # Test kinematic variables for different processes
        # Use s = 50.0 to ensure it's above threshold for all processes
        s = 50.0  # Mandelstam s [fm⁻²]
        M_squared = 1.0  # Arbitrary matrix element
        
        @testset "Identical quark scattering (uu → uu)" begin
            # u + u → u + u
            m1 = 1.52
            m2 = 1.52
            m3 = 1.52
            m4 = 1.52
            
            # Check threshold
            s_threshold = (m1 + m2)^2
            @test s > s_threshold
            @test check_kinematic_threshold(s, m1, m2, warn_close=false) == true
            
            # Calculate kinematic variables
            s_12_plus = s - (m1 + m2)^2
            s_12_minus = s - (m1 - m2)^2
            
            # Since m1 = m2, s_12_minus = s
            @test s_12_minus ≈ s
            
            # Calculate cross-section
            dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
            @test isfinite(dsigma_dt)
            @test dsigma_dt > 0.0
            
            # Verify formula
            expected = M_squared / (16π * s_12_plus * s_12_minus)
            @test dsigma_dt ≈ expected rtol=1e-12
        end
        
        @testset "Different quark scattering (ud → ud)" begin
            # u + d → u + d (m_u = m_d in our model)
            m1 = 1.52
            m2 = 1.52
            m3 = 1.52
            m4 = 1.52
            
            s_12_plus = s - (m1 + m2)^2
            s_12_minus = s - (m1 - m2)^2
            
            dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
            @test isfinite(dsigma_dt)
            @test dsigma_dt > 0.0
        end
        
        @testset "Strange quark scattering (us → us)" begin
            # u + s → u + s
            m1 = 1.52  # u quark
            m2 = 3.04  # s quark
            m3 = 1.52
            m4 = 3.04
            
            s_threshold = (m1 + m2)^2
            @test s > s_threshold
            @test check_kinematic_threshold(s, m1, m2, warn_close=false) == true
            
            s_12_plus = s - (m1 + m2)^2
            s_12_minus = s - (m1 - m2)^2
            
            # Since m1 ≠ m2, s_12_minus ≠ s
            @test s_12_minus < s
            @test s_12_minus > 0.0
            
            dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
            @test isfinite(dsigma_dt)
            @test dsigma_dt > 0.0
            
            expected = M_squared / (16π * s_12_plus * s_12_minus)
            @test dsigma_dt ≈ expected rtol=1e-12
        end
        
        @testset "Heavy quark scattering (ss → ss)" begin
            # s + s → s + s
            m1 = 3.04
            m2 = 3.04
            m3 = 3.04
            m4 = 3.04
            
            s_threshold = (m1 + m2)^2
            @test s > s_threshold
            @test check_kinematic_threshold(s, m1, m2, warn_close=false) == true
            
            s_12_plus = s - (m1 + m2)^2
            s_12_minus = s - (m1 - m2)^2
            
            # Since m1 = m2, s_12_minus = s
            @test s_12_minus ≈ s
            
            dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
            @test isfinite(dsigma_dt)
            @test dsigma_dt > 0.0
        end
        
        @testset "Quark-antiquark scattering (uubar → uubar)" begin
            # u + ubar → u + ubar
            m1 = 1.52
            m2 = 1.52  # Antiquark has same mass
            m3 = 1.52
            m4 = 1.52
            
            s_12_plus = s - (m1 + m2)^2
            s_12_minus = s - (m1 - m2)^2
            
            dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
            @test isfinite(dsigma_dt)
            @test dsigma_dt > 0.0
        end
        
        @testset "Flavor-changing process (uubar → ssbar)" begin
            # u + ubar → s + sbar
            m1 = 1.52  # u
            m2 = 1.52  # ubar
            m3 = 3.04  # s
            m4 = 3.04  # sbar
            
            # Initial state threshold
            s_threshold_initial = (m1 + m2)^2
            # Final state threshold
            s_threshold_final = (m3 + m4)^2
            
            # Need s > max(s_threshold_initial, s_threshold_final)
            s_high = max(s_threshold_initial, s_threshold_final) + 10.0
            
            @test check_kinematic_threshold(s_high, m1, m2, warn_close=false) == true
            @test check_kinematic_threshold(s_high, m3, m4, warn_close=false) == true
            
            s_12_plus = s_high - (m1 + m2)^2
            s_12_minus = s_high - (m1 - m2)^2
            
            dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
            @test isfinite(dsigma_dt)
            @test dsigma_dt > 0.0
        end
        
        @testset "Comparison between processes" begin
            # Compare cross-sections for different processes at same s
            s = 50.0
            M_squared = 1.0
            
            # uu → uu
            m_uu = 1.52
            s_12_plus_uu = s - (m_uu + m_uu)^2
            s_12_minus_uu = s - (m_uu - m_uu)^2
            dsigma_uu = differential_cross_section(s_12_plus_uu, s_12_minus_uu, M_squared)
            
            # us → us
            m_u = 1.52
            m_s = 3.04
            s_12_plus_us = s - (m_u + m_s)^2
            s_12_minus_us = s - (m_u - m_s)^2
            dsigma_us = differential_cross_section(s_12_plus_us, s_12_minus_us, M_squared)
            
            # ss → ss
            s_12_plus_ss = s - (m_s + m_s)^2
            s_12_minus_ss = s - (m_s - m_s)^2
            dsigma_ss = differential_cross_section(s_12_plus_ss, s_12_minus_ss, M_squared)
            
            # All should be positive and finite
            @test dsigma_uu > 0.0 && isfinite(dsigma_uu)
            @test dsigma_us > 0.0 && isfinite(dsigma_us)
            @test dsigma_ss > 0.0 && isfinite(dsigma_ss)
            
            # Heavier quarks have smaller s_12_plus (closer to threshold)
            # so larger cross-section (for same M²)
            @test s_12_plus_uu > s_12_plus_us > s_12_plus_ss
            @test dsigma_uu < dsigma_us < dsigma_ss
        end
        
        println("✓ Different scattering processes tests passed")
    end
    
    # ========================================================================
    # Test 4: Degenerate Cases and Numerical Stability
    # ========================================================================
    @testset "Degenerate cases and numerical stability" begin
        println("\n" * "="^70)
        println("Testing degenerate cases and numerical stability")
        println("="^70)
        
        @testset "Very small s_12_minus (m1 ≈ m2)" begin
            # When m1 ≈ m2, s_12_minus ≈ 0, which can cause numerical issues
            m1 = 1.52
            m2 = 1.52
            s = 20.0
            
            s_12_plus = s - (m1 + m2)^2
            s_12_minus = s - (m1 - m2)^2  # = s when m1 = m2
            M_squared = 1.0
            
            # Should handle gracefully
            dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
            @test isfinite(dsigma_dt)
            @test dsigma_dt > 0.0
        end
        
        @testset "Very small M_squared" begin
            # Test with very small matrix element
            m1 = 1.52
            m2 = 1.52
            s = 20.0
            
            s_12_plus = s - (m1 + m2)^2
            s_12_minus = s - (m1 - m2)^2
            
            for M_squared in [1e-10, 1e-8, 1e-6, 1e-4]
                dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
                @test isfinite(dsigma_dt)
                @test dsigma_dt >= 0.0
                
                # Verify formula
                expected = M_squared / (16π * s_12_plus * s_12_minus)
                @test dsigma_dt ≈ expected rtol=1e-12
            end
        end
        
        @testset "Zero M_squared" begin
            # When M² = 0, cross-section should be exactly 0
            m1 = 1.52
            m2 = 1.52
            s = 20.0
            
            s_12_plus = s - (m1 + m2)^2
            s_12_minus = s - (m1 - m2)^2
            M_squared = 0.0
            
            dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
            @test dsigma_dt == 0.0
        end
        
        @testset "Very large M_squared" begin
            # Test with very large matrix element
            m1 = 1.52
            m2 = 1.52
            s = 20.0
            
            s_12_plus = s - (m1 + m2)^2
            s_12_minus = s - (m1 - m2)^2
            
            for M_squared in [1e4, 1e6, 1e8, 1e10]
                dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
                @test isfinite(dsigma_dt)
                @test dsigma_dt > 0.0
                
                # Verify formula
                expected = M_squared / (16π * s_12_plus * s_12_minus)
                @test dsigma_dt ≈ expected rtol=1e-12
            end
        end
        
        @testset "Extreme mass ratios" begin
            # Test with very different masses
            m_light = 0.1
            m_heavy = 10.0
            s = 150.0  # Well above threshold
            
            s_threshold = (m_light + m_heavy)^2
            @test s > s_threshold
            
            s_12_plus = s - (m_light + m_heavy)^2
            s_12_minus = s - (m_light - m_heavy)^2
            M_squared = 1.0
            
            dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
            @test isfinite(dsigma_dt)
            @test dsigma_dt > 0.0
            
            # Verify formula
            expected = M_squared / (16π * s_12_plus * s_12_minus)
            @test dsigma_dt ≈ expected rtol=1e-12
        end
        
        println("✓ Degenerate cases and numerical stability tests passed")
    end
    
end

println("\n" * "="^70)
println("DifferentialCrossSection edge case tests complete!")
println("="^70)
println("\nSummary:")
println("  ✓ Threshold behavior tested (s near kinematic threshold)")
println("  ✓ High-energy behavior tested")
println("  ✓ Different scattering processes tested")
println("  ✓ Degenerate cases and numerical stability tested")
println("="^70)
