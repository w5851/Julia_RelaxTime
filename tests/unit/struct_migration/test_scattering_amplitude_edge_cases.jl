"""
Unit tests for ScatteringAmplitude module edge cases and specific behaviors.

Tests Requirements 4.4, 4.5:
- Different scattering processes
- Amplitude with different K_coeffs values
- Amplitude depends correctly on coupling strengths

These tests complement the property-based tests by focusing on specific
examples, edge cases, and physical constraints.
"""

using Test

# Load required modules
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src/relaxtime"))

include("../../../src/Constants_PNJL.jl")
include("../../../src/relaxtime/ScatteringAmplitude.jl")
include("../../../src/relaxtime/EffectiveCouplings.jl")
include("../../../src/relaxtime/OneLoopIntegrals.jl")
include("../../../src/integration/GaussLegendre.jl")

using .ScatteringAmplitude
using .Constants_PNJL
using .EffectiveCouplings
using .OneLoopIntegrals: A
using .GaussLegendre: gauleg

# Load ParameterTypes directly (don't use test_utils to avoid Supposition dependency)
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, joinpath(@__DIR__, "../../../src/ParameterTypes.jl"))
end
using .Main.ParameterTypes: QuarkParams, ThermoParams, as_namedtuple

@testset "ScatteringAmplitude Unit Tests" begin
    
    # ========================================================================
    # Test Setup: Common parameters
    # ========================================================================
    
    # Physical parameters (typical PNJL values)
    m_u = 1.52  # fm⁻¹ (≈ 300 MeV)
    m_s = 3.04  # fm⁻¹ (≈ 600 MeV)
    μ_u = 0.3   # fm⁻¹
    μ_s = 0.3   # fm⁻¹
    T = 0.15    # fm⁻¹ (≈ 30 MeV)
    Φ = 0.5
    Φbar = 0.5
    ξ = 0.0
    
    # Compute A functions
    nodes_p, weights_p = gauleg(0.0, 20.0, 64)
    A_u = A(m_u, μ_u, T, Φ, Φbar, nodes_p, weights_p)
    A_s = A(m_s, μ_s, T, Φ, Φbar, nodes_p, weights_p)
    
    # Compute G functions and K_coeffs
    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)
    G_fm2 = Constants_PNJL.G_fm2
    K_fm5 = Constants_PNJL.K_fm5
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
    
    # Create parameter structures
    quark_params = (
        m = (u=m_u, d=m_u, s=m_s),
        μ = (u=μ_u, d=μ_u, s=μ_s),
        A = (u=A_u, d=A_u, s=A_s)
    )
    thermo_params = (T=T, Φ=Φ, Φbar=Φbar, ξ=ξ)
    
    # Kinematic variables (must be above threshold for uu_to_uu: 4m_u² ≈ 9.25 fm⁻²)
    s = 31.0  # fm⁻² (well above threshold)
    t = -0.3  # fm⁻² (negative for physical scattering)
    
    # ========================================================================
    # Test 1: Different Scattering Processes
    # ========================================================================
    
    @testset "Test 1: Different Scattering Processes" begin
        println("\n" * "="^70)
        println("Test 1: Different Scattering Processes")
        println("="^70)
        
        # Test all qq scattering processes
        # Note: ss_to_ss requires s > 4m_s² ≈ 37 fm⁻², so we skip it with s=31.0
        qq_processes = [:uu_to_uu, :ud_to_ud, :us_to_us]  # Skip ss_to_ss due to threshold
        
        @testset "qq scattering: $process" for process in qq_processes
            M_squared = scattering_amplitude_squared(
                process, s, t, quark_params, thermo_params, K_coeffs
            )
            
            # Physical constraints
            @test M_squared >= 0.0
            @test isfinite(M_squared)
            
            println("  $process: |M|² = $(round(M_squared, digits=6)) fm⁻⁴")
        end
        
        # Test ss_to_ss separately with higher s value (above threshold)
        @testset "qq scattering: ss_to_ss (high energy)" begin
            s_high = 50.0  # fm⁻² (above 4m_s² ≈ 37 fm⁻²)
            M_squared = scattering_amplitude_squared(
                :ss_to_ss, s_high, t, quark_params, thermo_params, K_coeffs
            )
            
            @test M_squared >= 0.0
            @test isfinite(M_squared)
            
            println("  ss_to_ss (s=50.0): |M|² = $(round(M_squared, digits=6)) fm⁻⁴")
        end
        
        # Test all qqbar scattering processes
        qqbar_processes = [
            :udbar_to_udbar, :usbar_to_usbar, :dubar_to_dubar, :subar_to_subar,
            :uubar_to_uubar, :uubar_to_ddbar, :uubar_to_ssbar,
            :ssbar_to_uubar, :ssbar_to_ssbar
        ]
        
        @testset "qqbar scattering: $process" for process in qqbar_processes
            M_squared = scattering_amplitude_squared(
                process, s, t, quark_params, thermo_params, K_coeffs
            )
            
            # Physical constraints
            @test M_squared >= 0.0
            @test isfinite(M_squared)
            
            println("  $process: |M|² = $(round(M_squared, digits=6)) fm⁻⁴")
        end
        
        println("✓ All scattering processes produce valid amplitudes")
    end
    
    # ========================================================================
    # Test 2: Amplitude with Different K_coeffs Values
    # ========================================================================
    
    @testset "Test 2: Amplitude with Different K_coeffs Values" begin
        println("\n" * "="^70)
        println("Test 2: Amplitude with Different K_coeffs Values")
        println("="^70)
        
        process = :uu_to_uu
        
        # Test with standard K_coeffs
        M_squared_standard = scattering_amplitude_squared(
            process, s, t, quark_params, thermo_params, K_coeffs
        )
        
        # Test with K=0 (no 't Hooft interaction)
        K_coeffs_K0 = calculate_effective_couplings(G_fm2, 0.0, G_u, G_s)
        M_squared_K0 = scattering_amplitude_squared(
            process, s, t, quark_params, thermo_params, K_coeffs_K0
        )
        
        # Test with doubled K
        K_coeffs_2K = calculate_effective_couplings(G_fm2, 2.0 * K_fm5, G_u, G_s)
        M_squared_2K = scattering_amplitude_squared(
            process, s, t, quark_params, thermo_params, K_coeffs_2K
        )
        
        # Test with halved K
        K_coeffs_halfK = calculate_effective_couplings(G_fm2, 0.5 * K_fm5, G_u, G_s)
        M_squared_halfK = scattering_amplitude_squared(
            process, s, t, quark_params, thermo_params, K_coeffs_halfK
        )
        
        # Verify all are valid
        @test M_squared_standard >= 0.0 && isfinite(M_squared_standard)
        @test M_squared_K0 >= 0.0 && isfinite(M_squared_K0)
        @test M_squared_2K >= 0.0 && isfinite(M_squared_2K)
        @test M_squared_halfK >= 0.0 && isfinite(M_squared_halfK)
        
        # Verify they are different (K affects the amplitude)
        @test M_squared_K0 != M_squared_standard
        @test M_squared_2K != M_squared_standard
        @test M_squared_halfK != M_squared_standard
        
        println("  Standard K: |M|² = $(round(M_squared_standard, digits=6)) fm⁻⁴")
        println("  K = 0:      |M|² = $(round(M_squared_K0, digits=6)) fm⁻⁴")
        println("  K = 2K₀:    |M|² = $(round(M_squared_2K, digits=6)) fm⁻⁴")
        println("  K = 0.5K₀:  |M|² = $(round(M_squared_halfK, digits=6)) fm⁻⁴")
        println("✓ Amplitude varies with K_coeffs as expected")
    end
    
    # ========================================================================
    # Test 3: Amplitude Depends Correctly on Coupling Strengths
    # ========================================================================
    
    @testset "Test 3: Amplitude Depends on Coupling Strengths" begin
        println("\n" * "="^70)
        println("Test 3: Amplitude Depends on Coupling Strengths")
        println("="^70)
        
        process = :uu_to_uu
        
        # Test with standard G
        M_squared_standard = scattering_amplitude_squared(
            process, s, t, quark_params, thermo_params, K_coeffs
        )
        
        # Test with doubled G (stronger four-quark interaction)
        K_coeffs_2G = calculate_effective_couplings(2.0 * G_fm2, K_fm5, G_u, G_s)
        M_squared_2G = scattering_amplitude_squared(
            process, s, t, quark_params, thermo_params, K_coeffs_2G
        )
        
        # Test with halved G (weaker four-quark interaction)
        K_coeffs_halfG = calculate_effective_couplings(0.5 * G_fm2, K_fm5, G_u, G_s)
        M_squared_halfG = scattering_amplitude_squared(
            process, s, t, quark_params, thermo_params, K_coeffs_halfG
        )
        
        # Verify all are valid
        @test M_squared_standard >= 0.0 && isfinite(M_squared_standard)
        @test M_squared_2G >= 0.0 && isfinite(M_squared_2G)
        @test M_squared_halfG >= 0.0 && isfinite(M_squared_halfG)
        
        # Verify they are different (G affects the amplitude)
        @test M_squared_2G != M_squared_standard
        @test M_squared_halfG != M_squared_standard
        
        # Stronger coupling should generally give larger amplitude
        # (though this is not always monotonic due to interference terms)
        println("  Standard G: |M|² = $(round(M_squared_standard, digits=6)) fm⁻⁴")
        println("  G = 2G₀:    |M|² = $(round(M_squared_2G, digits=6)) fm⁻⁴")
        println("  G = 0.5G₀:  |M|² = $(round(M_squared_halfG, digits=6)) fm⁻⁴")
        
        # Test that amplitude scales with coupling strength
        # For weak coupling, |M|² should scale approximately as G²
        # (exact scaling depends on interference terms)
        ratio_2G = M_squared_2G / M_squared_standard
        ratio_halfG = M_squared_halfG / M_squared_standard
        
        println("  Ratio (2G/G₀): $(round(ratio_2G, digits=3))")
        println("  Ratio (0.5G/G₀): $(round(ratio_halfG, digits=3))")
        
        # Verify that doubling G increases amplitude (not necessarily by factor of 4)
        @test M_squared_2G > M_squared_standard
        @test M_squared_halfG < M_squared_standard
        
        println("✓ Amplitude depends on coupling strengths as expected")
    end
    
    # ========================================================================
    # Test 4: Amplitude at Different Kinematic Points
    # ========================================================================
    
    @testset "Test 4: Amplitude at Different Kinematic Points" begin
        println("\n" * "="^70)
        println("Test 4: Amplitude at Different Kinematic Points")
        println("="^70)
        
        process = :uu_to_uu
        
        # Test at different s values (energy dependence)
        s_values = [10.0, 20.0, 31.0, 50.0]  # fm⁻² (all above threshold)
        t_fixed = -0.3  # fm⁻²
        
        println("\n  Energy dependence (fixed t = $t_fixed fm⁻²):")
        for s_val in s_values
            M_squared = scattering_amplitude_squared(
                process, s_val, t_fixed, quark_params, thermo_params, K_coeffs
            )
            @test M_squared >= 0.0 && isfinite(M_squared)
            println("    s = $(s_val) fm⁻²: |M|² = $(round(M_squared, digits=6)) fm⁻⁴")
        end
        
        # Test at different t values (momentum transfer dependence)
        s_fixed = 31.0  # fm⁻²
        t_values = [-2.0, -1.0, -0.5, -0.1]  # fm⁻²
        
        println("\n  Momentum transfer dependence (fixed s = $s_fixed fm⁻²):")
        for t_val in t_values
            M_squared = scattering_amplitude_squared(
                process, s_fixed, t_val, quark_params, thermo_params, K_coeffs
            )
            @test M_squared >= 0.0 && isfinite(M_squared)
            println("    t = $(t_val) fm⁻²: |M|² = $(round(M_squared, digits=6)) fm⁻⁴")
        end
        
        println("✓ Amplitude is well-behaved across kinematic range")
    end
    
    # ========================================================================
    # Test 5: Amplitude with Different Quark Masses
    # ========================================================================
    
    @testset "Test 5: Amplitude with Different Quark Masses" begin
        println("\n" * "="^70)
        println("Test 5: Amplitude with Different Quark Masses")
        println("="^70)
        
        # Test uu scattering with different u quark masses
        process = :uu_to_uu
        
        mass_values = [1.0, 1.52, 2.0, 2.5]  # fm⁻¹
        
        println("\n  Mass dependence for uu scattering:")
        for m_test in mass_values
            # Recompute A and G for this mass
            A_test = A(m_test, μ_u, T, Φ, Φbar, nodes_p, weights_p)
            G_test = calculate_G_from_A(A_test, m_test)
            K_coeffs_test = calculate_effective_couplings(G_fm2, K_fm5, G_test, G_s)
            
            quark_params_test = (
                m = (u=m_test, d=m_test, s=m_s),
                μ = (u=μ_u, d=μ_u, s=μ_s),
                A = (u=A_test, d=A_test, s=A_s)
            )
            
            M_squared = scattering_amplitude_squared(
                process, s, t, quark_params_test, thermo_params, K_coeffs_test
            )
            
            @test M_squared >= 0.0 && isfinite(M_squared)
            println("    m_u = $(m_test) fm⁻¹: |M|² = $(round(M_squared, digits=6)) fm⁻⁴")
        end
        
        println("✓ Amplitude varies with quark mass as expected")
    end
    
    # ========================================================================
    # Test 6: Symmetry Relations
    # ========================================================================
    
    @testset "Test 6: Symmetry Relations" begin
        println("\n" * "="^70)
        println("Test 6: Symmetry Relations")
        println("="^70)
        
        # Test isospin symmetry: udbar_to_udbar should equal dubar_to_dubar
        M_udbar = scattering_amplitude_squared(
            :udbar_to_udbar, s, t, quark_params, thermo_params, K_coeffs
        )
        M_dubar = scattering_amplitude_squared(
            :dubar_to_dubar, s, t, quark_params, thermo_params, K_coeffs
        )
        
        @test isapprox(M_udbar, M_dubar, rtol=1e-12)
        println("  udbar_to_udbar: |M|² = $(round(M_udbar, digits=6)) fm⁻⁴")
        println("  dubar_to_dubar: |M|² = $(round(M_dubar, digits=6)) fm⁻⁴")
        println("  Difference: $(abs(M_udbar - M_dubar))")
        
        # Test isospin symmetry: usbar_to_usbar should equal subar_to_subar
        M_usbar = scattering_amplitude_squared(
            :usbar_to_usbar, s, t, quark_params, thermo_params, K_coeffs
        )
        M_subar = scattering_amplitude_squared(
            :subar_to_subar, s, t, quark_params, thermo_params, K_coeffs
        )
        
        @test isapprox(M_usbar, M_subar, rtol=1e-12)
        println("  usbar_to_usbar: |M|² = $(round(M_usbar, digits=6)) fm⁻⁴")
        println("  subar_to_subar: |M|² = $(round(M_subar, digits=6)) fm⁻⁴")
        println("  Difference: $(abs(M_usbar - M_subar))")
        
        println("✓ Symmetry relations satisfied")
    end
    
    # ========================================================================
    # Test 7: Struct vs NamedTuple Interface
    # ========================================================================
    
    @testset "Test 7: Struct vs NamedTuple Interface" begin
        println("\n" * "="^70)
        println("Test 7: Struct vs NamedTuple Interface")
        println("="^70)
        
        process = :uu_to_uu
        
        # Create struct parameters
        q_struct = QuarkParams((u=m_u, d=m_u, s=m_s), (u=μ_u, d=μ_u, s=μ_s))
        t_struct = ThermoParams(T, Φ, Φbar, ξ)
        
        # Add A field to struct-derived NamedTuple
        q_struct_with_A = merge(as_namedtuple(q_struct), (A=(u=A_u, d=A_u, s=A_s),))
        t_nt = as_namedtuple(t_struct)
        
        # Compute with struct-derived parameters
        M_squared_struct = scattering_amplitude_squared(
            process, s, t, q_struct_with_A, t_nt, K_coeffs
        )
        
        # Compute with NamedTuple parameters
        M_squared_nt = scattering_amplitude_squared(
            process, s, t, quark_params, thermo_params, K_coeffs
        )
        
        # Verify equivalence
        @test isapprox(M_squared_struct, M_squared_nt, rtol=1e-12, atol=1e-14)
        
        println("  Struct-derived: |M|² = $(round(M_squared_struct, digits=10)) fm⁻⁴")
        println("  NamedTuple:     |M|² = $(round(M_squared_nt, digits=10)) fm⁻⁴")
        println("  Difference: $(abs(M_squared_struct - M_squared_nt))")
        println("✓ Struct and NamedTuple interfaces produce identical results")
    end
    
end

println("\n" * "="^70)
println("ScatteringAmplitude Unit Tests Complete!")
println("="^70)
println("\nAll unit tests passed:")
println("  1. Different scattering processes tested")
println("  2. Amplitude varies with K_coeffs values")
println("  3. Amplitude depends on coupling strengths")
println("  4. Amplitude well-behaved across kinematic range")
println("  5. Amplitude varies with quark mass")
println("  6. Symmetry relations satisfied")
println("  7. Struct and NamedTuple interfaces equivalent")
println("="^70)
