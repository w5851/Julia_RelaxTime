"""
Unit tests for ParticleSymbols module edge cases.

Tests error handling and edge cases for:
- Invalid particle symbols
- All particle types (:u, :d, :s, :ubar, :dbar, :sbar)
- All scattering processes
- Requirements: 7.1, 7.2, 7.3
"""

using Test

# Load ParameterTypes module
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, joinpath(@__DIR__, "../../../src/ParameterTypes.jl"))
end

using Main.ParameterTypes: QuarkParams, ThermoParams, as_namedtuple

# Load ParticleSymbols module
if !isdefined(Main, :ParticleSymbols)
    Base.include(Main, joinpath(@__DIR__, "../../../src/utils/ParticleSymbols.jl"))
end

using Main.ParticleSymbols: get_mass, get_chemical_potential, get_quark_masses_for_process, extract_flavor

@testset "ParticleSymbols Edge Cases" begin
    
    # Create test parameters
    q_struct = QuarkParams((u=1.52, d=1.52, s=3.04), (u=0.3, d=0.3, s=0.3))
    q_nt = as_namedtuple(q_struct)
    
    # ========== Test Error Handling for Invalid Particle Symbols ==========
    @testset "Error handling for invalid particle symbols" begin
        # Test get_mass with invalid particle symbol
        @test_throws ErrorException get_mass(:invalid, q_struct)
        @test_throws ErrorException get_mass(:invalid, q_nt)
        
        # Test get_mass with unknown flavor (after extract_flavor)
        @test_throws ErrorException get_mass(:c, q_struct)  # charm quark not supported
        @test_throws ErrorException get_mass(:b, q_struct)  # bottom quark not supported
        
        # Test get_chemical_potential with invalid particle symbol
        @test_throws ErrorException get_chemical_potential(:invalid, q_struct)
        @test_throws ErrorException get_chemical_potential(:invalid, q_nt)
        
        # Test get_chemical_potential with unknown flavor
        @test_throws ErrorException get_chemical_potential(:c, q_struct)
        @test_throws ErrorException get_chemical_potential(:t, q_struct)
    end
    
    # ========== Test All Particle Types ==========
    @testset "All particle types with struct parameters" begin
        # Test quarks
        @test get_mass(:u, q_struct) ≈ 1.52
        @test get_mass(:d, q_struct) ≈ 1.52
        @test get_mass(:s, q_struct) ≈ 3.04
        
        # Test antiquarks (should have same mass as quarks)
        @test get_mass(:ubar, q_struct) ≈ 1.52
        @test get_mass(:dbar, q_struct) ≈ 1.52
        @test get_mass(:sbar, q_struct) ≈ 3.04
        
        # Test chemical potentials for quarks
        @test get_chemical_potential(:u, q_struct) ≈ 0.3
        @test get_chemical_potential(:d, q_struct) ≈ 0.3
        @test get_chemical_potential(:s, q_struct) ≈ 0.3
        
        # Test chemical potentials for antiquarks (should return same positive value)
        @test get_chemical_potential(:ubar, q_struct) ≈ 0.3
        @test get_chemical_potential(:dbar, q_struct) ≈ 0.3
        @test get_chemical_potential(:sbar, q_struct) ≈ 0.3
    end
    
    @testset "All particle types with NamedTuple parameters" begin
        # Test quarks
        @test get_mass(:u, q_nt) ≈ 1.52
        @test get_mass(:d, q_nt) ≈ 1.52
        @test get_mass(:s, q_nt) ≈ 3.04
        
        # Test antiquarks
        @test get_mass(:ubar, q_nt) ≈ 1.52
        @test get_mass(:dbar, q_nt) ≈ 1.52
        @test get_mass(:sbar, q_nt) ≈ 3.04
        
        # Test chemical potentials for quarks
        @test get_chemical_potential(:u, q_nt) ≈ 0.3
        @test get_chemical_potential(:d, q_nt) ≈ 0.3
        @test get_chemical_potential(:s, q_nt) ≈ 0.3
        
        # Test chemical potentials for antiquarks
        @test get_chemical_potential(:ubar, q_nt) ≈ 0.3
        @test get_chemical_potential(:dbar, q_nt) ≈ 0.3
        @test get_chemical_potential(:sbar, q_nt) ≈ 0.3
    end
    
    # ========== Test extract_flavor Function ==========
    @testset "extract_flavor edge cases" begin
        # Test quarks
        @test extract_flavor(:u) == (:u, false)
        @test extract_flavor(:d) == (:d, false)
        @test extract_flavor(:s) == (:s, false)
        
        # Test antiquarks
        @test extract_flavor(:ubar) == (:u, true)
        @test extract_flavor(:dbar) == (:d, true)
        @test extract_flavor(:sbar) == (:s, true)
        
        # Test that extract_flavor handles the "bar" suffix correctly
        # (This tests the fallback path in extract_flavor)
        flavor, is_anti = extract_flavor(:custombar)
        @test flavor == :custom
        @test is_anti == true
    end
    
    # ========== Test All Scattering Processes ==========
    @testset "All scattering processes with struct parameters" begin
        # Test quark-quark scattering processes
        @testset "Quark-quark scattering" begin
            # ud_to_ud: (m_u, m_d, m_u, m_d)
            masses = get_quark_masses_for_process(:ud_to_ud, q_struct)
            @test masses[1] ≈ 1.52  # u
            @test masses[2] ≈ 1.52  # d
            @test masses[3] ≈ 1.52  # u
            @test masses[4] ≈ 1.52  # d
            
            # uu_to_uu: (m_u, m_u, m_u, m_u)
            masses = get_quark_masses_for_process(:uu_to_uu, q_struct)
            @test all(m -> m ≈ 1.52, masses)
            
            # us_to_us: (m_u, m_s, m_u, m_s)
            masses = get_quark_masses_for_process(:us_to_us, q_struct)
            @test masses[1] ≈ 1.52  # u
            @test masses[2] ≈ 3.04  # s
            @test masses[3] ≈ 1.52  # u
            @test masses[4] ≈ 3.04  # s
            
            # ss_to_ss: (m_s, m_s, m_s, m_s)
            masses = get_quark_masses_for_process(:ss_to_ss, q_struct)
            @test all(m -> m ≈ 3.04, masses)
        end
        
        @testset "Antiquark-antiquark scattering" begin
            # ubardbar_to_ubardbar: (m_u, m_d, m_u, m_d)
            masses = get_quark_masses_for_process(:ubardbar_to_ubardbar, q_struct)
            @test masses[1] ≈ 1.52  # ubar
            @test masses[2] ≈ 1.52  # dbar
            @test masses[3] ≈ 1.52  # ubar
            @test masses[4] ≈ 1.52  # dbar
            
            # ubarubar_to_ubarubar: (m_u, m_u, m_u, m_u)
            masses = get_quark_masses_for_process(:ubarubar_to_ubarubar, q_struct)
            @test all(m -> m ≈ 1.52, masses)
            
            # ubarsbar_to_ubarsbar: (m_u, m_s, m_u, m_s)
            masses = get_quark_masses_for_process(:ubarsbar_to_ubarsbar, q_struct)
            @test masses[1] ≈ 1.52  # ubar
            @test masses[2] ≈ 3.04  # sbar
            @test masses[3] ≈ 1.52  # ubar
            @test masses[4] ≈ 3.04  # sbar
            
            # sbarsbar_to_sbarsbar: (m_s, m_s, m_s, m_s)
            masses = get_quark_masses_for_process(:sbarsbar_to_sbarsbar, q_struct)
            @test all(m -> m ≈ 3.04, masses)
        end
        
        @testset "Quark-antiquark scattering" begin
            # udbar_to_udbar: (m_u, m_d, m_u, m_d)
            masses = get_quark_masses_for_process(:udbar_to_udbar, q_struct)
            @test masses[1] ≈ 1.52  # u
            @test masses[2] ≈ 1.52  # dbar
            @test masses[3] ≈ 1.52  # u
            @test masses[4] ≈ 1.52  # dbar
            
            # dubar_to_dubar: (m_d, m_u, m_d, m_u)
            masses = get_quark_masses_for_process(:dubar_to_dubar, q_struct)
            @test masses[1] ≈ 1.52  # d
            @test masses[2] ≈ 1.52  # ubar
            @test masses[3] ≈ 1.52  # d
            @test masses[4] ≈ 1.52  # ubar
            
            # uubar_to_uubar: (m_u, m_u, m_u, m_u)
            masses = get_quark_masses_for_process(:uubar_to_uubar, q_struct)
            @test all(m -> m ≈ 1.52, masses)
            
            # uubar_to_ddbar: (m_u, m_u, m_d, m_d)
            masses = get_quark_masses_for_process(:uubar_to_ddbar, q_struct)
            @test masses[1] ≈ 1.52  # u
            @test masses[2] ≈ 1.52  # ubar
            @test masses[3] ≈ 1.52  # d
            @test masses[4] ≈ 1.52  # dbar
            
            # usbar_to_usbar: (m_u, m_s, m_u, m_s)
            masses = get_quark_masses_for_process(:usbar_to_usbar, q_struct)
            @test masses[1] ≈ 1.52  # u
            @test masses[2] ≈ 3.04  # sbar
            @test masses[3] ≈ 1.52  # u
            @test masses[4] ≈ 3.04  # sbar
            
            # subar_to_subar: (m_s, m_u, m_s, m_u)
            masses = get_quark_masses_for_process(:subar_to_subar, q_struct)
            @test masses[1] ≈ 3.04  # s
            @test masses[2] ≈ 1.52  # ubar
            @test masses[3] ≈ 3.04  # s
            @test masses[4] ≈ 1.52  # ubar
            
            # uubar_to_ssbar: (m_u, m_u, m_s, m_s)
            masses = get_quark_masses_for_process(:uubar_to_ssbar, q_struct)
            @test masses[1] ≈ 1.52  # u
            @test masses[2] ≈ 1.52  # ubar
            @test masses[3] ≈ 3.04  # s
            @test masses[4] ≈ 3.04  # sbar
            
            # ssbar_to_uubar: (m_s, m_s, m_u, m_u)
            masses = get_quark_masses_for_process(:ssbar_to_uubar, q_struct)
            @test masses[1] ≈ 3.04  # s
            @test masses[2] ≈ 3.04  # sbar
            @test masses[3] ≈ 1.52  # u
            @test masses[4] ≈ 1.52  # ubar
            
            # ssbar_to_ssbar: (m_s, m_s, m_s, m_s)
            masses = get_quark_masses_for_process(:ssbar_to_ssbar, q_struct)
            @test all(m -> m ≈ 3.04, masses)
        end
    end
    
    @testset "All scattering processes with NamedTuple parameters" begin
        # Test a representative sample with NamedTuple to ensure equivalence
        
        # Quark-quark
        masses_struct = get_quark_masses_for_process(:uu_to_uu, q_struct)
        masses_nt = get_quark_masses_for_process(:uu_to_uu, q_nt)
        @test all(i -> masses_struct[i] ≈ masses_nt[i], 1:4)
        
        # Antiquark-antiquark
        masses_struct = get_quark_masses_for_process(:ubarubar_to_ubarubar, q_struct)
        masses_nt = get_quark_masses_for_process(:ubarubar_to_ubarubar, q_nt)
        @test all(i -> masses_struct[i] ≈ masses_nt[i], 1:4)
        
        # Quark-antiquark
        masses_struct = get_quark_masses_for_process(:uubar_to_ssbar, q_struct)
        masses_nt = get_quark_masses_for_process(:uubar_to_ssbar, q_nt)
        @test all(i -> masses_struct[i] ≈ masses_nt[i], 1:4)
    end
    
    # ========== Test Edge Cases with Different Mass Values ==========
    @testset "Edge cases with different mass values" begin
        # Test with zero chemical potentials
        q_zero_mu = QuarkParams((u=1.52, d=1.52, s=3.04), (u=0.0, d=0.0, s=0.0))
        @test get_chemical_potential(:u, q_zero_mu) ≈ 0.0
        @test get_chemical_potential(:ubar, q_zero_mu) ≈ 0.0
        
        # Test with very small masses
        q_small = QuarkParams((u=0.1, d=0.1, s=0.1), (u=0.1, d=0.1, s=0.1))
        @test get_mass(:u, q_small) ≈ 0.1
        @test get_mass(:s, q_small) ≈ 0.1
        masses = get_quark_masses_for_process(:us_to_us, q_small)
        @test all(m -> m ≈ 0.1, masses)
        
        # Test with large masses
        q_large = QuarkParams((u=10.0, d=10.0, s=20.0), (u=1.0, d=1.0, s=1.0))
        @test get_mass(:s, q_large) ≈ 20.0
        masses = get_quark_masses_for_process(:uubar_to_ssbar, q_large)
        @test masses[1] ≈ 10.0  # u
        @test masses[2] ≈ 10.0  # ubar
        @test masses[3] ≈ 20.0  # s
        @test masses[4] ≈ 20.0  # sbar
        
        # Test with asymmetric masses (u ≠ d)
        q_asym = QuarkParams((u=1.0, d=2.0, s=3.0), (u=0.1, d=0.2, s=0.3))
        masses = get_quark_masses_for_process(:ud_to_ud, q_asym)
        @test masses[1] ≈ 1.0  # u
        @test masses[2] ≈ 2.0  # d
        @test masses[3] ≈ 1.0  # u
        @test masses[4] ≈ 2.0  # d
        
        # Verify chemical potentials are also correct
        @test get_chemical_potential(:u, q_asym) ≈ 0.1
        @test get_chemical_potential(:d, q_asym) ≈ 0.2
        @test get_chemical_potential(:s, q_asym) ≈ 0.3
    end
    
    # ========== Test Return Types and Properties ==========
    @testset "Return types and properties" begin
        # Test that get_mass returns Float64
        @test get_mass(:u, q_struct) isa Float64
        @test get_mass(:ubar, q_nt) isa Float64
        
        # Test that get_chemical_potential returns Float64
        @test get_chemical_potential(:u, q_struct) isa Float64
        @test get_chemical_potential(:sbar, q_nt) isa Float64
        
        # Test that get_quark_masses_for_process returns NTuple{4, Float64}
        masses = get_quark_masses_for_process(:uu_to_uu, q_struct)
        @test masses isa NTuple{4, Float64}
        @test length(masses) == 4
        
        # Test that all masses are positive
        @test all(>(0), masses)
        
        # Test that all masses are finite
        @test all(isfinite, masses)
    end
    
end

println("\n" * "="^70)
println("ParticleSymbols edge case tests complete!")
println("="^70)
