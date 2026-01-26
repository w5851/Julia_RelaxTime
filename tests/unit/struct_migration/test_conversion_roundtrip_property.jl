"""
Property-based test for conversion round-trip preservation.

Property 3: Conversion Round-Trip Preservation
Validates: Requirements 9.5

Tests that converting QuarkParams/ThermoParams to NamedTuple and back preserves all field values.
"""

using Test
using Supposition

# Load test utilities
include("test_utils.jl")

# Load ParameterTypes module
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, joinpath(@__DIR__, "../../../src/ParameterTypes.jl"))
end

using Main.ParameterTypes: QuarkParams, ThermoParams, as_namedtuple

@testset "Property: Conversion Round-Trip Preservation" begin
    
    @testset "QuarkParams round-trip preservation" begin
        # Feature: parameter-struct-migration, Property 3: Conversion Round-Trip Preservation
        # **Validates: Requirements 9.5**
        @check max_examples=20 function property_quark_params_roundtrip(
            m_u = Data.Floats{Float64}(minimum=0.1, maximum=5.0),
            m_d = Data.Floats{Float64}(minimum=0.1, maximum=5.0),
            m_s = Data.Floats{Float64}(minimum=0.1, maximum=10.0),
            μ_u = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
            μ_d = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
            μ_s = Data.Floats{Float64}(minimum=0.0, maximum=1.0)
        )
            # Create original struct
            q_original = QuarkParams(
                (u=m_u, d=m_d, s=m_s),  # m
                (u=μ_u, d=μ_d, s=μ_s)   # μ
            )
            
            # Convert to NamedTuple
            q_nt = as_namedtuple(q_original)
            
            # Reconstruct struct from NamedTuple
            q_reconstructed = QuarkParams(q_nt)
            
            # Verify all field values are preserved exactly
            q_reconstructed.m.u == q_original.m.u || error("m.u not preserved: $(q_reconstructed.m.u) != $(q_original.m.u)")
            q_reconstructed.m.d == q_original.m.d || error("m.d not preserved: $(q_reconstructed.m.d) != $(q_original.m.d)")
            q_reconstructed.m.s == q_original.m.s || error("m.s not preserved: $(q_reconstructed.m.s) != $(q_original.m.s)")
            q_reconstructed.μ.u == q_original.μ.u || error("μ.u not preserved: $(q_reconstructed.μ.u) != $(q_original.μ.u)")
            q_reconstructed.μ.d == q_original.μ.d || error("μ.d not preserved: $(q_reconstructed.μ.d) != $(q_original.μ.d)")
            q_reconstructed.μ.s == q_original.μ.s || error("μ.s not preserved: $(q_reconstructed.μ.s) != $(q_original.μ.s)")
            
            # Verify NamedTuple has correct structure
            haskey(q_nt, :m) || error("NamedTuple missing :m key")
            haskey(q_nt, :μ) || error("NamedTuple missing :μ key")
            haskey(q_nt.m, :u) || error("NamedTuple.m missing :u key")
            haskey(q_nt.m, :d) || error("NamedTuple.m missing :d key")
            haskey(q_nt.m, :s) || error("NamedTuple.m missing :s key")
            haskey(q_nt.μ, :u) || error("NamedTuple.μ missing :u key")
            haskey(q_nt.μ, :d) || error("NamedTuple.μ missing :d key")
            haskey(q_nt.μ, :s) || error("NamedTuple.μ missing :s key")
            
            # Verify NamedTuple values match original
            q_nt.m.u == q_original.m.u || error("NamedTuple m.u mismatch: $(q_nt.m.u) != $(q_original.m.u)")
            q_nt.m.d == q_original.m.d || error("NamedTuple m.d mismatch: $(q_nt.m.d) != $(q_original.m.d)")
            q_nt.m.s == q_original.m.s || error("NamedTuple m.s mismatch: $(q_nt.m.s) != $(q_original.m.s)")
            q_nt.μ.u == q_original.μ.u || error("NamedTuple μ.u mismatch: $(q_nt.μ.u) != $(q_original.μ.u)")
            q_nt.μ.d == q_original.μ.d || error("NamedTuple μ.d mismatch: $(q_nt.μ.d) != $(q_original.μ.d)")
            q_nt.μ.s == q_original.μ.s || error("NamedTuple μ.s mismatch: $(q_nt.μ.s) != $(q_original.μ.s)")
        end
    end
    
    @testset "ThermoParams round-trip preservation" begin
        # Feature: parameter-struct-migration, Property 3: Conversion Round-Trip Preservation
        # **Validates: Requirements 9.5**
        @check max_examples=20 function property_thermo_params_roundtrip(
            T = Data.Floats{Float64}(minimum=0.05, maximum=0.3),
            Φ = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
            Φbar = Data.Floats{Float64}(minimum=0.0, maximum=1.0)
        )
            # Create original struct (ξ defaults to 0.0)
            t_original = ThermoParams(T, Φ, Φbar, 0.0)
            
            # Convert to NamedTuple
            t_nt = as_namedtuple(t_original)
            
            # Reconstruct struct from NamedTuple
            t_reconstructed = ThermoParams(t_nt)
            
            # Verify all field values are preserved exactly
            t_reconstructed.T == t_original.T || error("T not preserved: $(t_reconstructed.T) != $(t_original.T)")
            t_reconstructed.Φ == t_original.Φ || error("Φ not preserved: $(t_reconstructed.Φ) != $(t_original.Φ)")
            t_reconstructed.Φbar == t_original.Φbar || error("Φbar not preserved: $(t_reconstructed.Φbar) != $(t_original.Φbar)")
            t_reconstructed.ξ == t_original.ξ || error("ξ not preserved: $(t_reconstructed.ξ) != $(t_original.ξ)")
            
            # Verify NamedTuple has correct structure
            haskey(t_nt, :T) || error("NamedTuple missing :T key")
            haskey(t_nt, :Φ) || error("NamedTuple missing :Φ key")
            haskey(t_nt, :Φbar) || error("NamedTuple missing :Φbar key")
            haskey(t_nt, :ξ) || error("NamedTuple missing :ξ key")
            
            # Verify NamedTuple values match original
            t_nt.T == t_original.T || error("NamedTuple T mismatch: $(t_nt.T) != $(t_original.T)")
            t_nt.Φ == t_original.Φ || error("NamedTuple Φ mismatch: $(t_nt.Φ) != $(t_original.Φ)")
            t_nt.Φbar == t_original.Φbar || error("NamedTuple Φbar mismatch: $(t_nt.Φbar) != $(t_original.Φbar)")
            t_nt.ξ == t_original.ξ || error("NamedTuple ξ mismatch: $(t_nt.ξ) != $(t_original.ξ)")
        end
    end
    
    @testset "ThermoParams with non-zero ξ round-trip" begin
        # Feature: parameter-struct-migration, Property 3: Conversion Round-Trip Preservation
        # **Validates: Requirements 9.5**
        @check max_examples=20 function property_thermo_params_with_xi_roundtrip(
            T = Data.Floats{Float64}(minimum=0.05, maximum=0.3),
            Φ = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
            Φbar = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
            ξ = Data.Floats{Float64}(minimum=-0.8, maximum=0.8)
        )
            # Create original struct with non-zero ξ
            t_original = ThermoParams(T, Φ, Φbar, ξ)
            
            # Convert to NamedTuple
            t_nt = as_namedtuple(t_original)
            
            # Reconstruct struct from NamedTuple
            t_reconstructed = ThermoParams(t_nt)
            
            # Verify all field values are preserved exactly
            t_reconstructed.T == t_original.T || error("T not preserved: $(t_reconstructed.T) != $(t_original.T)")
            t_reconstructed.Φ == t_original.Φ || error("Φ not preserved: $(t_reconstructed.Φ) != $(t_original.Φ)")
            t_reconstructed.Φbar == t_original.Φbar || error("Φbar not preserved: $(t_reconstructed.Φbar) != $(t_original.Φbar)")
            t_reconstructed.ξ == t_original.ξ || error("ξ not preserved: $(t_reconstructed.ξ) != $(t_original.ξ)")
            
            # Verify NamedTuple values match original
            t_nt.T == t_original.T || error("NamedTuple T mismatch: $(t_nt.T) != $(t_original.T)")
            t_nt.Φ == t_original.Φ || error("NamedTuple Φ mismatch: $(t_nt.Φ) != $(t_original.Φ)")
            t_nt.Φbar == t_original.Φbar || error("NamedTuple Φbar mismatch: $(t_nt.Φbar) != $(t_original.Φbar)")
            t_nt.ξ == t_original.ξ || error("NamedTuple ξ mismatch: $(t_nt.ξ) != $(t_original.ξ)")
        end
    end
    
    @testset "Multiple round-trips preserve values" begin
        # Test that multiple conversions don't accumulate errors
        # Feature: parameter-struct-migration, Property 3: Conversion Round-Trip Preservation
        # **Validates: Requirements 9.5**
        @check max_examples=10 function property_multiple_roundtrips(
            m_u = Data.Floats{Float64}(minimum=0.1, maximum=5.0),
            m_d = Data.Floats{Float64}(minimum=0.1, maximum=5.0),
            m_s = Data.Floats{Float64}(minimum=0.1, maximum=10.0),
            μ_u = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
            μ_d = Data.Floats{Float64}(minimum=0.0, maximum=1.0),
            μ_s = Data.Floats{Float64}(minimum=0.0, maximum=1.0)
        )
            # Create original struct
            q_original = QuarkParams(
                (u=m_u, d=m_d, s=m_s),  # m
                (u=μ_u, d=μ_d, s=μ_s)   # μ
            )
            
            # Perform multiple round-trips
            q_current = q_original
            for i in 1:5
                q_nt = as_namedtuple(q_current)
                q_current = QuarkParams(q_nt)
            end
            
            # Verify values are still preserved after multiple conversions
            q_current.m.u == q_original.m.u || error("m.u changed after multiple roundtrips: $(q_current.m.u) != $(q_original.m.u)")
            q_current.m.d == q_original.m.d || error("m.d changed after multiple roundtrips: $(q_current.m.d) != $(q_original.m.d)")
            q_current.m.s == q_original.m.s || error("m.s changed after multiple roundtrips: $(q_current.m.s) != $(q_original.m.s)")
            q_current.μ.u == q_original.μ.u || error("μ.u changed after multiple roundtrips: $(q_current.μ.u) != $(q_original.μ.u)")
            q_current.μ.d == q_original.μ.d || error("μ.d changed after multiple roundtrips: $(q_current.μ.d) != $(q_original.μ.d)")
            q_current.μ.s == q_original.μ.s || error("μ.s changed after multiple roundtrips: $(q_current.μ.s) != $(q_original.μ.s)")
        end
    end
    
    @testset "Edge cases: extreme values" begin
        # Test round-trip with extreme but valid values
        @testset "Very small masses" begin
            q = QuarkParams((u=0.01, d=0.01, s=0.01), (u=0.0, d=0.0, s=0.0))
            q_nt = as_namedtuple(q)
            q_reconstructed = QuarkParams(q_nt)
            
            @test q_reconstructed.m.u == q.m.u
            @test q_reconstructed.m.d == q.m.d
            @test q_reconstructed.m.s == q.m.s
        end
        
        @testset "Very large masses" begin
            q = QuarkParams((u=10.0, d=10.0, s=20.0), (u=0.5, d=0.5, s=0.5))
            q_nt = as_namedtuple(q)
            q_reconstructed = QuarkParams(q_nt)
            
            @test q_reconstructed.m.u == q.m.u
            @test q_reconstructed.m.d == q.m.d
            @test q_reconstructed.m.s == q.m.s
        end
        
        @testset "Zero chemical potentials" begin
            q = QuarkParams((u=1.5, d=1.5, s=3.0), (u=0.0, d=0.0, s=0.0))
            q_nt = as_namedtuple(q)
            q_reconstructed = QuarkParams(q_nt)
            
            @test q_reconstructed.μ.u == 0.0
            @test q_reconstructed.μ.d == 0.0
            @test q_reconstructed.μ.s == 0.0
        end
        
        @testset "Very low temperature" begin
            t = ThermoParams(0.01, 0.1, 0.1, 0.0)
            t_nt = as_namedtuple(t)
            t_reconstructed = ThermoParams(t_nt)
            
            @test t_reconstructed.T == t.T
        end
        
        @testset "Very high temperature" begin
            t = ThermoParams(0.5, 0.9, 0.9, 0.0)
            t_nt = as_namedtuple(t)
            t_reconstructed = ThermoParams(t_nt)
            
            @test t_reconstructed.T == t.T
        end
    end
end
