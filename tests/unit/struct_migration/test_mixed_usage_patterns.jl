# Test mixed struct/NamedTuple usage patterns
# Validates Requirement 8.3: Mixed struct/NamedTuple usage

using Test

# Load test utilities
include("test_utils.jl")

# Load modules
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "..", "src", "ParameterTypes.jl"))
end

using Main.ParameterTypes: QuarkParams, ThermoParams, as_namedtuple

# Load ParticleSymbols module first (required by other modules)
if !isdefined(Main, :ParticleSymbols)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "..", "src", "utils", "ParticleSymbols.jl"))
end

# Load ScatteringAmplitude module
if !isdefined(Main, :ScatteringAmplitude)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "ScatteringAmplitude.jl"))
end

# Load DifferentialCrossSection module
if !isdefined(Main, :DifferentialCrossSection)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "DifferentialCrossSection.jl"))
end

# Load TotalPropagator module
if !isdefined(Main, :TotalPropagator)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "TotalPropagator.jl"))
end

# Load TotalCrossSection module
if !isdefined(Main, :TotalCrossSection)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "TotalCrossSection.jl"))
end

# Load AverageScatteringRate module
if !isdefined(Main, :AverageScatteringRate)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "AverageScatteringRate.jl"))
end

# Load RelaxationTime module
if !isdefined(Main, :RelaxationTime)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxationTime.jl"))
end

@testset "Mixed Struct/NamedTuple Usage Patterns" begin
    
    # Create test parameters in both formats
    q_struct = QuarkParams(
        (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
    )
    t_struct = ThermoParams((T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0))
    
    q_nt = as_namedtuple(q_struct)
    t_nt = as_namedtuple(t_struct)
    
    K_coeffs = (K_pi=1.0, K_K=1.0, K_eta=1.0, K_sigma=1.0)
    
    @testset "Pattern 1: Struct quark_params + NamedTuple thermo_params" begin
        # Test with TotalCrossSection
        s = 5.0
        xs = Main.TotalCrossSection.total_cross_section(
            :us_to_us,
            s,
            q_struct,  # Struct
            t_nt,      # NamedTuple
            K_coeffs
        )
        
        @test isfinite(xs)
        @test xs >= 0
        
        # Test with ParticleSymbols
        m_u = Main.ParticleSymbols.get_mass(:u, q_struct)  # Struct
        @test m_u == 1.52
        
        mu_u = Main.ParticleSymbols.get_chemical_potential(:u, q_struct)  # Struct
        @test mu_u == 0.3
    end
    
    @testset "Pattern 2: NamedTuple quark_params + Struct thermo_params" begin
        # Test with TotalCrossSection
        s = 5.0
        xs = Main.TotalCrossSection.total_cross_section(
            :us_to_us,
            s,
            q_nt,      # NamedTuple
            t_struct,  # Struct
            K_coeffs
        )
        
        @test isfinite(xs)
        @test xs >= 0
        
        # Test with ParticleSymbols
        m_u = Main.ParticleSymbols.get_mass(:u, q_nt)  # NamedTuple
        @test m_u == 1.52
        
        mu_u = Main.ParticleSymbols.get_chemical_potential(:u, q_nt)  # NamedTuple
        @test mu_u == 0.3
    end
    
    @testset "Pattern 3: Mixed usage in workflow" begin
        # Simulate a workflow that mixes both representations
        s = 5.0
        
        # Step 1: Use struct parameters
        xs1 = Main.TotalCrossSection.total_cross_section(
            :us_to_us,
            s,
            q_struct,
            t_struct,
            K_coeffs
        )
        
        # Step 2: Convert to NamedTuple and use
        q_nt_converted = as_namedtuple(q_struct)
        t_nt_converted = as_namedtuple(t_struct)
        
        xs2 = Main.TotalCrossSection.total_cross_section(
            :us_to_us,
            s,
            q_nt_converted,
            t_nt_converted,
            K_coeffs
        )
        
        # Step 3: Mix struct and NamedTuple
        xs3 = Main.TotalCrossSection.total_cross_section(
            :us_to_us,
            s,
            q_struct,
            t_nt_converted,
            K_coeffs
        )
        
        # All results should be identical
        @test isapprox(xs1, xs2, rtol=1e-12)
        @test isapprox(xs1, xs3, rtol=1e-12)
    end
    
    @testset "Pattern 4: Verify all combinations work correctly" begin
        # Test all four combinations for a simple function
        s = 5.0
        
        # Struct + Struct
        xs_ss = Main.TotalCrossSection.total_cross_section(
            :us_to_us, s, q_struct, t_struct, K_coeffs
        )
        
        # Struct + NamedTuple
        xs_sn = Main.TotalCrossSection.total_cross_section(
            :us_to_us, s, q_struct, t_nt, K_coeffs
        )
        
        # NamedTuple + Struct
        xs_ns = Main.TotalCrossSection.total_cross_section(
            :us_to_us, s, q_nt, t_struct, K_coeffs
        )
        
        # NamedTuple + NamedTuple
        xs_nn = Main.TotalCrossSection.total_cross_section(
            :us_to_us, s, q_nt, t_nt, K_coeffs
        )
        
        # All should be identical
        @test isapprox(xs_ss, xs_sn, rtol=1e-12)
        @test isapprox(xs_ss, xs_ns, rtol=1e-12)
        @test isapprox(xs_ss, xs_nn, rtol=1e-12)
        @test isapprox(xs_sn, xs_ns, rtol=1e-12)
        @test isapprox(xs_sn, xs_nn, rtol=1e-12)
        @test isapprox(xs_ns, xs_nn, rtol=1e-12)
    end
    
    @testset "Pattern 5: Nested function calls with mixed types" begin
        # Test that mixed types work correctly through nested function calls
        s = 5.0
        
        # Call total_cross_section (which calls scattering_amplitude internally)
        # with struct quark_params and NamedTuple thermo_params
        xs1 = Main.TotalCrossSection.total_cross_section(
            :us_to_us, s, q_struct, t_nt, K_coeffs
        )
        
        # Call with NamedTuple quark_params and struct thermo_params
        xs2 = Main.TotalCrossSection.total_cross_section(
            :us_to_us, s, q_nt, t_struct, K_coeffs
        )
        
        # Both should be identical
        @test isapprox(xs1, xs2, rtol=1e-12)
    end
end

println("\n" * "="^80)
println("Mixed Struct/NamedTuple Usage Pattern Tests Complete!")
println("="^80)
println("✓ Pattern 1: Struct quark_params + NamedTuple thermo_params")
println("✓ Pattern 2: NamedTuple quark_params + Struct thermo_params")
println("✓ Pattern 3: Mixed usage in workflow")
println("✓ Pattern 4: All combinations work correctly")
println("✓ Pattern 5: Nested function calls with mixed types")
println("="^80)
