"""
Property-based tests for ParticleSymbols module struct support.

Tests Property 2: Field Extraction Correctness
**Validates: Requirements 7.1, 7.2, 7.3**

Tests:
1. get_mass with random QuarkParams structs vs NamedTuples
2. get_chemical_potential with random QuarkParams structs vs NamedTuples
3. get_quark_masses_for_process with random parameters
4. Verify extracted values match struct fields exactly
"""

using Test
using Supposition

# Load test utilities
include("test_utils.jl")

# Load ParticleSymbols module
if !isdefined(Main, :ParticleSymbols)
    Base.include(Main, joinpath(@__DIR__, "../../../src/utils/ParticleSymbols.jl"))
end

using Main.ParticleSymbols: get_mass, get_chemical_potential, get_quark_masses_for_process

@testset "ParticleSymbols Struct Support" begin
    
    # ========== Property 2: Field Extraction Correctness ==========
    @testset "Property 2: Field Extraction Correctness" begin
        # Feature: parameter-struct-migration, Property 2: Field Extraction Correctness
        # Validates: Requirements 7.1, 7.2, 7.3
        
        # Test get_mass for all particle types
        @testset "get_mass equivalence" begin
            @check max_examples=100 function property_get_mass_equivalence(
                m_u = Data.Floats{Float64}(minimum=0.1, maximum=5.0),
                m_d = Data.Floats{Float64}(minimum=0.1, maximum=5.0),
                m_s = Data.Floats{Float64}(minimum=0.1, maximum=10.0),
                μ_u = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
                μ_d = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
                μ_s = Data.Floats{Float64}(minimum=0.0, maximum=1.0)
            )
                # Create struct and NamedTuple parameters
                q_struct = QuarkParams((u=m_u, d=m_d, s=m_s), (u=μ_u, d=μ_d, s=μ_s))
                q_nt = as_namedtuple(q_struct)
                
                # Test u quark
                mass_struct_u = get_mass(:u, q_struct)
                mass_nt_u = get_mass(:u, q_nt)
                mass_struct_u ≈ mass_nt_u || error("get_mass(:u) not equivalent")
                mass_struct_u ≈ m_u || error("get_mass(:u) doesn't match m.u")
                
                # Test d quark
                mass_struct_d = get_mass(:d, q_struct)
                mass_nt_d = get_mass(:d, q_nt)
                mass_struct_d ≈ mass_nt_d || error("get_mass(:d) not equivalent")
                mass_struct_d ≈ m_d || error("get_mass(:d) doesn't match m.d")
                
                # Test s quark
                mass_struct_s = get_mass(:s, q_struct)
                mass_nt_s = get_mass(:s, q_nt)
                mass_struct_s ≈ mass_nt_s || error("get_mass(:s) not equivalent")
                mass_struct_s ≈ m_s || error("get_mass(:s) doesn't match m.s")
                
                # Test antiquarks (should have same mass as quarks)
                mass_struct_ubar = get_mass(:ubar, q_struct)
                mass_nt_ubar = get_mass(:ubar, q_nt)
                mass_struct_ubar ≈ mass_nt_ubar || error("get_mass(:ubar) not equivalent")
                mass_struct_ubar ≈ m_u || error("get_mass(:ubar) doesn't match m.u")
                
                mass_struct_dbar = get_mass(:dbar, q_struct)
                mass_nt_dbar = get_mass(:dbar, q_nt)
                mass_struct_dbar ≈ mass_nt_dbar || error("get_mass(:dbar) not equivalent")
                mass_struct_dbar ≈ m_d || error("get_mass(:dbar) doesn't match m.d")
                
                mass_struct_sbar = get_mass(:sbar, q_struct)
                mass_nt_sbar = get_mass(:sbar, q_nt)
                mass_struct_sbar ≈ mass_nt_sbar || error("get_mass(:sbar) not equivalent")
                mass_struct_sbar ≈ m_s || error("get_mass(:sbar) doesn't match m.s")
            end
        end
        
        # Test get_chemical_potential for all particle types
        @testset "get_chemical_potential equivalence" begin
            @check max_examples=100 function property_get_chemical_potential_equivalence(
                m_u = Data.Floats{Float64}(minimum=0.1, maximum=5.0),
                m_d = Data.Floats{Float64}(minimum=0.1, maximum=5.0),
                m_s = Data.Floats{Float64}(minimum=0.1, maximum=10.0),
                μ_u = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
                μ_d = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
                μ_s = Data.Floats{Float64}(minimum=0.0, maximum=1.0)
            )
                # Create struct and NamedTuple parameters
                q_struct = QuarkParams((u=m_u, d=m_d, s=m_s), (u=μ_u, d=μ_d, s=μ_s))
                q_nt = as_namedtuple(q_struct)
                
                # Test u quark
                μ_struct_u = get_chemical_potential(:u, q_struct)
                μ_nt_u = get_chemical_potential(:u, q_nt)
                μ_struct_u ≈ μ_nt_u || error("get_chemical_potential(:u) not equivalent")
                μ_struct_u ≈ μ_u || error("get_chemical_potential(:u) doesn't match μ.u")
                
                # Test d quark
                μ_struct_d = get_chemical_potential(:d, q_struct)
                μ_nt_d = get_chemical_potential(:d, q_nt)
                μ_struct_d ≈ μ_nt_d || error("get_chemical_potential(:d) not equivalent")
                μ_struct_d ≈ μ_d || error("get_chemical_potential(:d) doesn't match μ.d")
                
                # Test s quark
                μ_struct_s = get_chemical_potential(:s, q_struct)
                μ_nt_s = get_chemical_potential(:s, q_nt)
                μ_struct_s ≈ μ_nt_s || error("get_chemical_potential(:s) not equivalent")
                μ_struct_s ≈ μ_s || error("get_chemical_potential(:s) doesn't match μ.s")
                
                # Test antiquarks (should return same positive chemical potential)
                μ_struct_ubar = get_chemical_potential(:ubar, q_struct)
                μ_nt_ubar = get_chemical_potential(:ubar, q_nt)
                μ_struct_ubar ≈ μ_nt_ubar || error("get_chemical_potential(:ubar) not equivalent")
                μ_struct_ubar ≈ μ_u || error("get_chemical_potential(:ubar) doesn't match μ.u")
                
                μ_struct_dbar = get_chemical_potential(:dbar, q_struct)
                μ_nt_dbar = get_chemical_potential(:dbar, q_nt)
                μ_struct_dbar ≈ μ_nt_dbar || error("get_chemical_potential(:dbar) not equivalent")
                μ_struct_dbar ≈ μ_d || error("get_chemical_potential(:dbar) doesn't match μ.d")
                
                μ_struct_sbar = get_chemical_potential(:sbar, q_struct)
                μ_nt_sbar = get_chemical_potential(:sbar, q_nt)
                μ_struct_sbar ≈ μ_nt_sbar || error("get_chemical_potential(:sbar) not equivalent")
                μ_struct_sbar ≈ μ_s || error("get_chemical_potential(:sbar) doesn't match μ.s")
            end
        end
        
        # Test get_quark_masses_for_process for various scattering processes
        @testset "get_quark_masses_for_process equivalence" begin
            @check max_examples=100 function property_get_quark_masses_for_process_equivalence(
                m_u = Data.Floats{Float64}(minimum=0.1, maximum=5.0),
                m_d = Data.Floats{Float64}(minimum=0.1, maximum=5.0),
                m_s = Data.Floats{Float64}(minimum=0.1, maximum=10.0),
                μ_u = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
                μ_d = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
                μ_s = Data.Floats{Float64}(minimum=0.0, maximum=1.0)
            )
                # Create struct and NamedTuple parameters
                q_struct = QuarkParams((u=m_u, d=m_d, s=m_s), (u=μ_u, d=μ_d, s=μ_s))
                q_nt = as_namedtuple(q_struct)
                
                # Test uu_to_uu process
                masses_struct = get_quark_masses_for_process(:uu_to_uu, q_struct)
                masses_nt = get_quark_masses_for_process(:uu_to_uu, q_nt)
                masses_struct[1] ≈ masses_nt[1] || error("uu_to_uu[1] not equivalent")
                masses_struct[2] ≈ masses_nt[2] || error("uu_to_uu[2] not equivalent")
                masses_struct[3] ≈ masses_nt[3] || error("uu_to_uu[3] not equivalent")
                masses_struct[4] ≈ masses_nt[4] || error("uu_to_uu[4] not equivalent")
                all(isfinite, masses_struct) || error("uu_to_uu returned non-finite values")
                all(>(0), masses_struct) || error("uu_to_uu returned non-positive values")
                
                # Test ud_to_ud process
                masses_struct = get_quark_masses_for_process(:ud_to_ud, q_struct)
                masses_nt = get_quark_masses_for_process(:ud_to_ud, q_nt)
                masses_struct[1] ≈ masses_nt[1] || error("ud_to_ud[1] not equivalent")
                masses_struct[2] ≈ masses_nt[2] || error("ud_to_ud[2] not equivalent")
                masses_struct[3] ≈ masses_nt[3] || error("ud_to_ud[3] not equivalent")
                masses_struct[4] ≈ masses_nt[4] || error("ud_to_ud[4] not equivalent")
                
                # Test uubar_to_ddbar process
                masses_struct = get_quark_masses_for_process(:uubar_to_ddbar, q_struct)
                masses_nt = get_quark_masses_for_process(:uubar_to_ddbar, q_nt)
                masses_struct[1] ≈ masses_nt[1] || error("uubar_to_ddbar[1] not equivalent")
                masses_struct[2] ≈ masses_nt[2] || error("uubar_to_ddbar[2] not equivalent")
                masses_struct[3] ≈ masses_nt[3] || error("uubar_to_ddbar[3] not equivalent")
                masses_struct[4] ≈ masses_nt[4] || error("uubar_to_ddbar[4] not equivalent")
                
                # Test uubar_to_ssbar process
                masses_struct = get_quark_masses_for_process(:uubar_to_ssbar, q_struct)
                masses_nt = get_quark_masses_for_process(:uubar_to_ssbar, q_nt)
                masses_struct[1] ≈ masses_nt[1] || error("uubar_to_ssbar[1] not equivalent")
                masses_struct[2] ≈ masses_nt[2] || error("uubar_to_ssbar[2] not equivalent")
                masses_struct[3] ≈ masses_nt[3] || error("uubar_to_ssbar[3] not equivalent")
                masses_struct[4] ≈ masses_nt[4] || error("uubar_to_ssbar[4] not equivalent")
            end
        end
        
        # Test that get_quark_masses_for_process returns correct masses for specific processes
        @testset "get_quark_masses_for_process correctness" begin
            @check max_examples=100 function property_get_quark_masses_for_process_correctness(
                m_u = Data.Floats{Float64}(minimum=0.1, maximum=5.0),
                m_d = Data.Floats{Float64}(minimum=0.1, maximum=5.0),
                m_s = Data.Floats{Float64}(minimum=0.1, maximum=10.0),
                μ_u = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
                μ_d = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
                μ_s = Data.Floats{Float64}(minimum=0.0, maximum=1.0)
            )
                # Create struct parameters
                q_struct = QuarkParams((u=m_u, d=m_d, s=m_s), (u=μ_u, d=μ_d, s=μ_s))
                
                # Test specific processes and verify correct masses are returned
                # Process: uu_to_uu should return (m_u, m_u, m_u, m_u)
                masses = get_quark_masses_for_process(:uu_to_uu, q_struct)
                masses[1] ≈ m_u || error("uu_to_uu: mass[1] should be m_u")
                masses[2] ≈ m_u || error("uu_to_uu: mass[2] should be m_u")
                masses[3] ≈ m_u || error("uu_to_uu: mass[3] should be m_u")
                masses[4] ≈ m_u || error("uu_to_uu: mass[4] should be m_u")
                
                # Process: ud_to_ud should return (m_u, m_d, m_u, m_d)
                masses = get_quark_masses_for_process(:ud_to_ud, q_struct)
                masses[1] ≈ m_u || error("ud_to_ud: mass[1] should be m_u")
                masses[2] ≈ m_d || error("ud_to_ud: mass[2] should be m_d")
                masses[3] ≈ m_u || error("ud_to_ud: mass[3] should be m_u")
                masses[4] ≈ m_d || error("ud_to_ud: mass[4] should be m_d")
                
                # Process: us_to_us should return (m_u, m_s, m_u, m_s)
                masses = get_quark_masses_for_process(:us_to_us, q_struct)
                masses[1] ≈ m_u || error("us_to_us: mass[1] should be m_u")
                masses[2] ≈ m_s || error("us_to_us: mass[2] should be m_s")
                masses[3] ≈ m_u || error("us_to_us: mass[3] should be m_u")
                masses[4] ≈ m_s || error("us_to_us: mass[4] should be m_s")
                
                # Process: uubar_to_ddbar should return (m_u, m_u, m_d, m_d)
                # Note: antiquarks have same mass as quarks
                masses = get_quark_masses_for_process(:uubar_to_ddbar, q_struct)
                masses[1] ≈ m_u || error("uubar_to_ddbar: mass[1] should be m_u")
                masses[2] ≈ m_u || error("uubar_to_ddbar: mass[2] should be m_u")
                masses[3] ≈ m_d || error("uubar_to_ddbar: mass[3] should be m_d")
                masses[4] ≈ m_d || error("uubar_to_ddbar: mass[4] should be m_d")
                
                # Process: uubar_to_ssbar should return (m_u, m_u, m_s, m_s)
                masses = get_quark_masses_for_process(:uubar_to_ssbar, q_struct)
                masses[1] ≈ m_u || error("uubar_to_ssbar: mass[1] should be m_u")
                masses[2] ≈ m_u || error("uubar_to_ssbar: mass[2] should be m_u")
                masses[3] ≈ m_s || error("uubar_to_ssbar: mass[3] should be m_s")
                masses[4] ≈ m_s || error("uubar_to_ssbar: mass[4] should be m_s")
            end
        end
    end
    
end

println("\n" * "="^70)
println("ParticleSymbols struct support tests complete!")
println("="^70)
