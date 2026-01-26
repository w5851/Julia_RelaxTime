"""
Tests for QuarkParams and ThermoParams struct construction and conversion.

Tests:
1. QuarkParams construction from NamedTuple
2. ThermoParams construction from NamedTuple
3. Struct to NamedTuple conversion (as_namedtuple)
4. Round-trip conversion preservation (Property 3)
5. Error handling for missing fields
6. Field access and validation
"""

using Test
using Supposition

# Load test utilities
include("test_utils.jl")

@testset "ParameterTypes Module" begin
    
    # ========== Test 1: QuarkParams Construction ==========
    @testset "QuarkParams construction" begin
        # Valid construction from NamedTuple
        q_nt = (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
        q = QuarkParams(q_nt)
        
        @test q.m.u ≈ 1.52
        @test q.m.d ≈ 1.52
        @test q.m.s ≈ 3.04
        @test q.μ.u ≈ 0.3
        @test q.μ.d ≈ 0.3
        @test q.μ.s ≈ 0.3
        
        # Direct construction
        q2 = QuarkParams((u=1.52, d=1.52, s=3.04), (u=0.3, d=0.3, s=0.3))
        @test q2.m.u ≈ 1.52
        @test q2.μ.u ≈ 0.3
    end
    
    # ========== Test 2: ThermoParams Construction ==========
    @testset "ThermoParams construction" begin
        # Valid construction from NamedTuple
        t_nt = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
        t = ThermoParams(t_nt)
        
        @test t.T ≈ 0.15
        @test t.Φ ≈ 0.5
        @test t.Φbar ≈ 0.5
        @test t.ξ ≈ 0.0
        
        # Direct construction
        t2 = ThermoParams(0.15, 0.5, 0.5, 0.0)
        @test t2.T ≈ 0.15
        @test t2.Φ ≈ 0.5
        
        # Construction with default ξ
        t_nt_no_xi = (T=0.15, Φ=0.5, Φbar=0.5)
        t3 = ThermoParams(t_nt_no_xi)
        @test t3.ξ ≈ 0.0  # Should default to 0.0
    end
    
    # ========== Test 3: Struct to NamedTuple Conversion ==========
    @testset "as_namedtuple conversion" begin
        # QuarkParams conversion
        q = QuarkParams((u=1.52, d=1.52, s=3.04), (u=0.3, d=0.3, s=0.3))
        q_nt = as_namedtuple(q)
        
        @test q_nt isa NamedTuple
        @test haskey(q_nt, :m)
        @test haskey(q_nt, :μ)
        @test q_nt.m.u ≈ 1.52
        @test q_nt.μ.u ≈ 0.3
        
        # ThermoParams conversion
        t = ThermoParams(0.15, 0.5, 0.5, 0.0)
        t_nt = as_namedtuple(t)
        
        @test t_nt isa NamedTuple
        @test haskey(t_nt, :T)
        @test haskey(t_nt, :Φ)
        @test haskey(t_nt, :Φbar)
        @test haskey(t_nt, :ξ)
        @test t_nt.T ≈ 0.15
        @test t_nt.Φ ≈ 0.5
    end
    
    # ========== Test 4: Property 3 - Round-Trip Conversion ==========
    @testset "Property 3: Round-trip conversion preservation" begin
        # Feature: parameter-struct-migration, Property 3: Conversion Round-Trip Preservation
        # Validates: Requirements 9.5
        
        @check max_examples=100 function property_quark_params_roundtrip(
            m_u = Data.Floats{Float64}(minimum=0.1, maximum=5.0),
            m_d = Data.Floats{Float64}(minimum=0.1, maximum=5.0),
            m_s = Data.Floats{Float64}(minimum=0.1, maximum=10.0),
            μ_u = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
            μ_d = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
            μ_s = Data.Floats{Float64}(minimum=0.0, maximum=1.0)
        )
            # Create original struct
            q_original = QuarkParams((u=m_u, d=m_d, s=m_s), (u=μ_u, d=μ_d, s=μ_s))
            
            # Convert to NamedTuple and back
            q_nt = as_namedtuple(q_original)
            q_reconstructed = QuarkParams(q_nt)
            
            # Verify all fields are preserved (using Supposition's assertion mechanism)
            q_reconstructed.m.u ≈ q_original.m.u || error("m.u not preserved")
            q_reconstructed.m.d ≈ q_original.m.d || error("m.d not preserved")
            q_reconstructed.m.s ≈ q_original.m.s || error("m.s not preserved")
            q_reconstructed.μ.u ≈ q_original.μ.u || error("μ.u not preserved")
            q_reconstructed.μ.d ≈ q_original.μ.d || error("μ.d not preserved")
            q_reconstructed.μ.s ≈ q_original.μ.s || error("μ.s not preserved")
        end
        
        @check max_examples=100 function property_thermo_params_roundtrip(
            T = Data.Floats{Float64}(minimum=0.05, maximum=0.3),
            Φ = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
            Φbar = Data.Floats{Float64}(minimum=0.0, maximum=1.0)
        )
            # Create original struct
            t_original = ThermoParams(T, Φ, Φbar, 0.0)
            
            # Convert to NamedTuple and back
            t_nt = as_namedtuple(t_original)
            t_reconstructed = ThermoParams(t_nt)
            
            # Verify all fields are preserved (using Supposition's assertion mechanism)
            t_reconstructed.T ≈ t_original.T || error("T not preserved")
            t_reconstructed.Φ ≈ t_original.Φ || error("Φ not preserved")
            t_reconstructed.Φbar ≈ t_original.Φbar || error("Φbar not preserved")
            t_reconstructed.ξ ≈ t_original.ξ || error("ξ not preserved")
        end
    end
    
    # ========== Test 5: Error Handling ==========
    @testset "Error handling for missing fields" begin
        # QuarkParams missing m field
        @test_throws ErrorException QuarkParams((μ=(u=0.3, d=0.3, s=0.3),))
        
        # QuarkParams missing μ field
        @test_throws ErrorException QuarkParams((m=(u=1.52, d=1.52, s=3.04),))
        
        # ThermoParams missing T field
        @test_throws ErrorException ThermoParams((Φ=0.5, Φbar=0.5, ξ=0.0))
        
        # ThermoParams missing Φ field
        @test_throws ErrorException ThermoParams((T=0.15, Φbar=0.5, ξ=0.0))
        
        # ThermoParams missing Φbar field
        @test_throws ErrorException ThermoParams((T=0.15, Φ=0.5, ξ=0.0))
    end
    
    # ========== Test 6: Field Access ==========
    @testset "Field access and validation" begin
        q = QuarkParams((u=1.52, d=1.52, s=3.04), (u=0.3, d=0.3, s=0.3))
        
        # Test that fields are accessible
        @test q.m isa NamedTuple
        @test q.μ isa NamedTuple
        @test keys(q.m) == (:u, :d, :s)
        @test keys(q.μ) == (:u, :d, :s)
        
        t = ThermoParams(0.15, 0.5, 0.5, 0.0)
        
        # Test that fields are accessible
        @test t.T isa Float64
        @test t.Φ isa Float64
        @test t.Φbar isa Float64
        @test t.ξ isa Float64
    end
    
end

println("\n" * "="^70)
println("ParameterTypes tests complete!")
println("="^70)
