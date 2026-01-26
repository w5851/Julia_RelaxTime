"""
Property-based test for ensure_quark_params_has_A function.

Property 5: Extended QuarkParams with A Field
Validates: Requirements 9.4

Tests that ensure_quark_params_has_A correctly handles QuarkParams with and without A field,
and works with both struct and NamedTuple inputs.
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

# Load RelaxationTime module to access ensure_quark_params_has_A
if !isdefined(Main, :RelaxationTime)
    Base.include(Main, joinpath(@__DIR__, "../../../src/relaxtime/RelaxationTime.jl"))
end

@testset "Property: Extended QuarkParams with A Field" begin
    
    @testset "QuarkParams without A field - struct input" begin
        # Feature: parameter-struct-migration, Property 5: Extended QuarkParams with A Field
        # **Validates: Requirements 9.4**
        @check max_examples=20 function property_struct_workflow_support(
            m_u = Data.Floats{Float64}(minimum=0.5, maximum=3.0),
            m_d = Data.Floats{Float64}(minimum=0.5, maximum=3.0),
            m_s = Data.Floats{Float64}(minimum=1.0, maximum=5.0),
            μ_u = Data.Floats{Float64}(minimum=0.0, maximum=0.5),
            μ_d = Data.Floats{Float64}(minimum=0.0, maximum=0.5),
            μ_s = Data.Floats{Float64}(minimum=0.0, maximum=0.5),
            T = Data.Floats{Float64}(minimum=0.1, maximum=0.25)
        )
            # Create struct parameters without A field
            q_struct = QuarkParams(
                (u=m_u, d=m_d, s=m_s),  # m
                (u=μ_u, d=μ_d, s=μ_s)   # μ
            )
            t_struct = ThermoParams(T, 0.5, 0.5, 0.0)
            
            # Call ensure_quark_params_has_A
            q_with_A = Main.RelaxationTime.ensure_quark_params_has_A(q_struct, t_struct)
            
            # Verify output is a NamedTuple
            q_with_A isa NamedTuple || error("Output is not a NamedTuple: $(typeof(q_with_A))")
            
            # Verify output contains all original fields
            haskey(q_with_A, :m) || error("Output missing :m key")
            haskey(q_with_A, :μ) || error("Output missing :μ key")
            haskey(q_with_A.m, :u) || error("Output.m missing :u key")
            haskey(q_with_A.m, :d) || error("Output.m missing :d key")
            haskey(q_with_A.m, :s) || error("Output.m missing :s key")
            haskey(q_with_A.μ, :u) || error("Output.μ missing :u key")
            haskey(q_with_A.μ, :d) || error("Output.μ missing :d key")
            haskey(q_with_A.μ, :s) || error("Output.μ missing :s key")
            
            # Verify original values are preserved
            q_with_A.m.u == m_u || error("m.u not preserved: $(q_with_A.m.u) != $m_u")
            q_with_A.m.d == m_d || error("m.d not preserved: $(q_with_A.m.d) != $m_d")
            q_with_A.m.s == m_s || error("m.s not preserved: $(q_with_A.m.s) != $m_s")
            q_with_A.μ.u == μ_u || error("μ.u not preserved: $(q_with_A.μ.u) != $μ_u")
            q_with_A.μ.d == μ_d || error("μ.d not preserved: $(q_with_A.μ.d) != $μ_d")
            q_with_A.μ.s == μ_s || error("μ.s not preserved: $(q_with_A.μ.s) != $μ_s")
            
            # Verify A field was added
            haskey(q_with_A, :A) || error("Output missing :A key")
            haskey(q_with_A.A, :u) || error("Output.A missing :u key")
            haskey(q_with_A.A, :d) || error("Output.A missing :d key")
            haskey(q_with_A.A, :s) || error("Output.A missing :s key")
            
            # Verify A values are finite
            isfinite(q_with_A.A.u) || error("A.u is not finite: $(q_with_A.A.u)")
            isfinite(q_with_A.A.d) || error("A.d is not finite: $(q_with_A.A.d)")
            isfinite(q_with_A.A.s) || error("A.s is not finite: $(q_with_A.A.s)")
        end
    end
    
    @testset "QuarkParams without A field - NamedTuple input" begin
        # Feature: parameter-struct-migration, Property 5: Extended QuarkParams with A Field
        # **Validates: Requirements 9.4**
        @check max_examples=20 function property_namedtuple_workflow_support(
            m_u = Data.Floats{Float64}(minimum=0.5, maximum=3.0),
            m_d = Data.Floats{Float64}(minimum=0.5, maximum=3.0),
            m_s = Data.Floats{Float64}(minimum=1.0, maximum=5.0),
            μ_u = Data.Floats{Float64}(minimum=0.0, maximum=0.5),
            μ_d = Data.Floats{Float64}(minimum=0.0, maximum=0.5),
            μ_s = Data.Floats{Float64}(minimum=0.0, maximum=0.5),
            T = Data.Floats{Float64}(minimum=0.1, maximum=0.25)
        )
            # Create NamedTuple parameters without A field
            q_nt = (
                m = (u=m_u, d=m_d, s=m_s),
                μ = (u=μ_u, d=μ_d, s=μ_s)
            )
            t_nt = (T=T, Φ=0.5, Φbar=0.5, ξ=0.0)
            
            # Call ensure_quark_params_has_A
            q_with_A = Main.RelaxationTime.ensure_quark_params_has_A(q_nt, t_nt)
            
            # Verify output is a NamedTuple
            q_with_A isa NamedTuple || error("Output is not a NamedTuple: $(typeof(q_with_A))")
            
            # Verify output contains all original fields
            haskey(q_with_A, :m) || error("Output missing :m key")
            haskey(q_with_A, :μ) || error("Output missing :μ key")
            
            # Verify original values are preserved
            q_with_A.m.u == m_u || error("m.u not preserved: $(q_with_A.m.u) != $m_u")
            q_with_A.m.d == m_d || error("m.d not preserved: $(q_with_A.m.d) != $m_d")
            q_with_A.m.s == m_s || error("m.s not preserved: $(q_with_A.m.s) != $m_s")
            q_with_A.μ.u == μ_u || error("μ.u not preserved: $(q_with_A.μ.u) != $μ_u")
            q_with_A.μ.d == μ_d || error("μ.d not preserved: $(q_with_A.μ.d) != $μ_d")
            q_with_A.μ.s == μ_s || error("μ.s not preserved: $(q_with_A.μ.s) != $μ_s")
            
            # Verify A field was added
            haskey(q_with_A, :A) || error("Output missing :A key")
            haskey(q_with_A.A, :u) || error("Output.A missing :u key")
            haskey(q_with_A.A, :d) || error("Output.A missing :d key")
            haskey(q_with_A.A, :s) || error("Output.A missing :s key")
            
            # Verify A values are finite
            isfinite(q_with_A.A.u) || error("A.u is not finite: $(q_with_A.A.u)")
            isfinite(q_with_A.A.d) || error("A.d is not finite: $(q_with_A.A.d)")
            isfinite(q_with_A.A.s) || error("A.s is not finite: $(q_with_A.A.s)")
        end
    end
    
    @testset "QuarkParams with A field - should preserve it" begin
        # Test that if A field already exists, it is preserved
        @testset "Struct input with A field" begin
            q_struct = QuarkParams(
                (u=1.5, d=1.5, s=3.0),  # m
                (u=0.3, d=0.3, s=0.3)   # μ
            )
            t_struct = ThermoParams(0.15, 0.5, 0.5, 0.0)
            
            # First call to add A field
            q_with_A = Main.RelaxationTime.ensure_quark_params_has_A(q_struct, t_struct)
            
            # Store original A values
            A_u_original = q_with_A.A.u
            A_d_original = q_with_A.A.d
            A_s_original = q_with_A.A.s
            
            # Second call with q_with_A (which already has A field)
            q_with_A_again = Main.RelaxationTime.ensure_quark_params_has_A(q_with_A, t_struct)
            
            # Verify A field was preserved (not recomputed)
            @test q_with_A_again.A.u == A_u_original
            @test q_with_A_again.A.d == A_d_original
            @test q_with_A_again.A.s == A_s_original
            
            # Verify other fields are still present
            @test q_with_A_again.m.u == 1.5
            @test q_with_A_again.m.d == 1.5
            @test q_with_A_again.m.s == 3.0
            @test q_with_A_again.μ.u == 0.3
            @test q_with_A_again.μ.d == 0.3
            @test q_with_A_again.μ.s == 0.3
        end
        
        @testset "NamedTuple input with A field" begin
            q_nt = (
                m = (u=1.5, d=1.5, s=3.0),
                μ = (u=0.3, d=0.3, s=0.3),
                A = (u=0.1, d=0.1, s=0.2)  # Pre-existing A field
            )
            t_nt = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
            
            # Call ensure_quark_params_has_A
            q_result = Main.RelaxationTime.ensure_quark_params_has_A(q_nt, t_nt)
            
            # Verify A field was preserved
            @test q_result.A.u == 0.1
            @test q_result.A.d == 0.1
            @test q_result.A.s == 0.2
            
            # Verify other fields are still present
            @test q_result.m.u == 1.5
            @test q_result.m.d == 1.5
            @test q_result.m.s == 3.0
        end
    end
    
    @testset "Struct and NamedTuple equivalence" begin
        # Feature: parameter-struct-migration, Property 5: Extended QuarkParams with A Field
        # **Validates: Requirements 9.4**
        @check max_examples=10 function property_struct_namedtuple_equivalence(
            m_u = Data.Floats{Float64}(minimum=0.5, maximum=3.0),
            m_d = Data.Floats{Float64}(minimum=0.5, maximum=3.0),
            m_s = Data.Floats{Float64}(minimum=1.0, maximum=5.0),
            μ_u = Data.Floats{Float64}(minimum=0.0, maximum=0.5),
            μ_d = Data.Floats{Float64}(minimum=0.0, maximum=0.5),
            μ_s = Data.Floats{Float64}(minimum=0.0, maximum=0.5),
            T = Data.Floats{Float64}(minimum=0.1, maximum=0.25)
        )
            # Create struct parameters
            q_struct = QuarkParams(
                (u=m_u, d=m_d, s=m_s),  # m
                (u=μ_u, d=μ_d, s=μ_s)   # μ
            )
            t_struct = ThermoParams(T, 0.5, 0.5, 0.0)
            
            # Create equivalent NamedTuple parameters
            q_nt = as_namedtuple(q_struct)
            t_nt = as_namedtuple(t_struct)
            
            # Call with struct parameters
            q_struct_result = Main.RelaxationTime.ensure_quark_params_has_A(q_struct, t_struct)
            
            # Call with NamedTuple parameters
            q_nt_result = Main.RelaxationTime.ensure_quark_params_has_A(q_nt, t_nt)
            
            # Verify both produce identical results
            q_struct_result.m.u == q_nt_result.m.u || error("m.u mismatch: $(q_struct_result.m.u) != $(q_nt_result.m.u)")
            q_struct_result.m.d == q_nt_result.m.d || error("m.d mismatch: $(q_struct_result.m.d) != $(q_nt_result.m.d)")
            q_struct_result.m.s == q_nt_result.m.s || error("m.s mismatch: $(q_struct_result.m.s) != $(q_nt_result.m.s)")
            q_struct_result.μ.u == q_nt_result.μ.u || error("μ.u mismatch: $(q_struct_result.μ.u) != $(q_nt_result.μ.u)")
            q_struct_result.μ.d == q_nt_result.μ.d || error("μ.d mismatch: $(q_struct_result.μ.d) != $(q_nt_result.μ.d)")
            q_struct_result.μ.s == q_nt_result.μ.s || error("μ.s mismatch: $(q_struct_result.μ.s) != $(q_nt_result.μ.s)")
            
            # Verify A values are identical (within floating-point tolerance)
            approx_equal(q_struct_result.A.u, q_nt_result.A.u) || error("A.u not approximately equal: $(q_struct_result.A.u) vs $(q_nt_result.A.u)")
            approx_equal(q_struct_result.A.d, q_nt_result.A.d) || error("A.d not approximately equal: $(q_struct_result.A.d) vs $(q_nt_result.A.d)")
            approx_equal(q_struct_result.A.s, q_nt_result.A.s) || error("A.s not approximately equal: $(q_struct_result.A.s) vs $(q_nt_result.A.s)")
        end
    end
    
    @testset "Edge cases" begin
        @testset "Very small masses" begin
            q = QuarkParams((u=0.1, d=0.1, s=0.2), (u=0.0, d=0.0, s=0.0))
            t = ThermoParams(0.15, 0.5, 0.5, 0.0)
            
            q_with_A = Main.RelaxationTime.ensure_quark_params_has_A(q, t)
            
            @test haskey(q_with_A, :A)
            @test isfinite(q_with_A.A.u)
            @test isfinite(q_with_A.A.d)
            @test isfinite(q_with_A.A.s)
        end
        
        @testset "Large masses" begin
            q = QuarkParams((u=5.0, d=5.0, s=10.0), (u=0.5, d=0.5, s=0.5))
            t = ThermoParams(0.15, 0.5, 0.5, 0.0)
            
            q_with_A = Main.RelaxationTime.ensure_quark_params_has_A(q, t)
            
            @test haskey(q_with_A, :A)
            @test isfinite(q_with_A.A.u)
            @test isfinite(q_with_A.A.d)
            @test isfinite(q_with_A.A.s)
        end
        
        @testset "Low temperature" begin
            q = QuarkParams((u=1.5, d=1.5, s=3.0), (u=0.2, d=0.2, s=0.2))
            t = ThermoParams(0.05, 0.3, 0.3, 0.0)
            
            q_with_A = Main.RelaxationTime.ensure_quark_params_has_A(q, t)
            
            @test haskey(q_with_A, :A)
            @test isfinite(q_with_A.A.u)
            @test isfinite(q_with_A.A.d)
            @test isfinite(q_with_A.A.s)
        end
        
        @testset "High temperature" begin
            q = QuarkParams((u=1.5, d=1.5, s=3.0), (u=0.4, d=0.4, s=0.4))
            t = ThermoParams(0.25, 0.7, 0.7, 0.0)
            
            q_with_A = Main.RelaxationTime.ensure_quark_params_has_A(q, t)
            
            @test haskey(q_with_A, :A)
            @test isfinite(q_with_A.A.u)
            @test isfinite(q_with_A.A.d)
            @test isfinite(q_with_A.A.s)
        end
        
        @testset "Non-zero ξ" begin
            q = QuarkParams((u=1.5, d=1.5, s=3.0), (u=0.3, d=0.3, s=0.3))
            t = ThermoParams(0.15, 0.5, 0.5, 0.2)
            
            q_with_A = Main.RelaxationTime.ensure_quark_params_has_A(q, t)
            
            @test haskey(q_with_A, :A)
            @test isfinite(q_with_A.A.u)
            @test isfinite(q_with_A.A.d)
            @test isfinite(q_with_A.A.s)
        end
    end
    
    @testset "Output structure validation" begin
        q = QuarkParams((u=1.5, d=1.5, s=3.0), (u=0.3, d=0.3, s=0.3))
        t = ThermoParams(0.15, 0.5, 0.5, 0.0)
        
        q_with_A = Main.RelaxationTime.ensure_quark_params_has_A(q, t)
        
        # Verify output is a NamedTuple (not a struct)
        @test q_with_A isa NamedTuple
        
        # Verify it has exactly 3 top-level keys: m, μ, A
        @test length(keys(q_with_A)) == 3
        @test :m in keys(q_with_A)
        @test :μ in keys(q_with_A)
        @test :A in keys(q_with_A)
        
        # Verify nested structure
        @test q_with_A.m isa NamedTuple
        @test q_with_A.μ isa NamedTuple
        @test q_with_A.A isa NamedTuple
        
        # Verify each nested tuple has u, d, s keys
        @test Set(keys(q_with_A.m)) == Set([:u, :d, :s])
        @test Set(keys(q_with_A.μ)) == Set([:u, :d, :s])
        @test Set(keys(q_with_A.A)) == Set([:u, :d, :s])
        
        # Verify all values are Float64
        @test q_with_A.m.u isa Float64
        @test q_with_A.μ.u isa Float64
        @test q_with_A.A.u isa Float64
    end
end
