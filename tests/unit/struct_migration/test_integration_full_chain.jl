"""
Integration tests for full call chain with struct parameters.

Tests that struct parameters flow correctly through the entire RelaxationTime module chain:
RelaxationTime → AverageScatteringRate → TotalCrossSection → ScatteringAmplitude → DifferentialCrossSection

Validates Requirements 5.4, 8.3
"""

using Test

# Load test utilities
include("test_utils.jl")

# Load source modules
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, joinpath(@__DIR__, "../../../src/ParameterTypes.jl"))
end

using Main.ParameterTypes: QuarkParams, ThermoParams, as_namedtuple

# Load RelaxationTime module (which includes all dependencies)
if !isdefined(Main, :RelaxationTime)
    Base.include(Main, joinpath(@__DIR__, "../../../src/relaxtime/RelaxationTime.jl"))
end

@testset "Integration: Full Call Chain" begin
    
    @testset "Struct parameters flow through entire chain" begin
        # Create struct parameters
        q_struct = QuarkParams(
            (u=1.52, d=1.52, s=3.04),  # m
            (u=0.3, d=0.3, s=0.3)       # μ
        )
        t_struct = ThermoParams(0.15, 0.5, 0.5, 0.0)
        
        # Create K_coeffs (coupling coefficients)
        K_coeffs = (
            K_pi = 1.0,
            K_K = 1.0,
            K_eta = 1.0,
            K_sigma = 1.0,
            K_kappa = 1.0,
            K_a0 = 1.0
        )
        
        # Create test densities
        densities = (
            u = 0.1,
            d = 0.1,
            s = 0.05,
            ubar = 0.1,
            dbar = 0.1,
            sbar = 0.05
        )
        
        # Call top-level function with struct parameters
        result = Main.RelaxationTime.relaxation_times(
            q_struct, t_struct, K_coeffs;
            densities = densities
        )
        
        # Verify result structure
        @test haskey(result, :tau)
        @test haskey(result, :tau_inv)
        @test haskey(result, :rates)
        
        # Verify tau has expected flavors
        @test haskey(result.tau, :u)
        @test haskey(result.tau, :d)
        @test haskey(result.tau, :s)
        
        # Verify all relaxation times are finite and positive
        @test isfinite(result.tau.u)
        @test isfinite(result.tau.d)
        @test isfinite(result.tau.s)
        @test result.tau.u > 0
        @test result.tau.d > 0
        @test result.tau.s > 0
        
        # Verify inverse relaxation times
        @test isfinite(result.tau_inv.u)
        @test isfinite(result.tau_inv.d)
        @test isfinite(result.tau_inv.s)
        @test result.tau_inv.u > 0
        @test result.tau_inv.d > 0
        @test result.tau_inv.s > 0
        
        # Verify rates structure
        @test result.rates isa NamedTuple
        @test length(result.rates) > 0
    end
    
    @testset "Struct vs NamedTuple equivalence through full chain" begin
        # Create struct parameters
        q_struct = QuarkParams(
            (u=1.52, d=1.52, s=3.04),  # m
            (u=0.3, d=0.3, s=0.3)       # μ
        )
        t_struct = ThermoParams(0.15, 0.5, 0.5, 0.0)
        
        # Create equivalent NamedTuple parameters
        q_nt = as_namedtuple(q_struct)
        t_nt = as_namedtuple(t_struct)
        
        K_coeffs = (
            K_pi = 1.0,
            K_K = 1.0,
            K_eta = 1.0,
            K_sigma = 1.0,
            K_kappa = 1.0,
            K_a0 = 1.0
        )
        
        densities = (
            u = 0.1,
            d = 0.1,
            s = 0.05,
            ubar = 0.1,
            dbar = 0.1,
            sbar = 0.05
        )
        
        # Call with struct parameters
        result_struct = Main.RelaxationTime.relaxation_times(
            q_struct, t_struct, K_coeffs;
            densities = densities
        )
        
        # Call with NamedTuple parameters
        result_nt = Main.RelaxationTime.relaxation_times(
            q_nt, t_nt, K_coeffs;
            densities = densities
        )
        
        # Verify equivalence of relaxation times
        @test approx_equal(result_struct.tau.u, result_nt.tau.u)
        @test approx_equal(result_struct.tau.d, result_nt.tau.d)
        @test approx_equal(result_struct.tau.s, result_nt.tau.s)
        
        # Verify equivalence of inverse relaxation times
        @test approx_equal(result_struct.tau_inv.u, result_nt.tau_inv.u)
        @test approx_equal(result_struct.tau_inv.d, result_nt.tau_inv.d)
        @test approx_equal(result_struct.tau_inv.s, result_nt.tau_inv.s)
        
        # Verify equivalence of rates
        @test keys(result_struct.rates) == keys(result_nt.rates)
        for process in keys(result_struct.rates)
            @test approx_equal(result_struct.rates[process], result_nt.rates[process])
        end
    end
    
    @testset "Different parameter combinations" begin
        K_coeffs = (
            K_pi = 1.0,
            K_K = 1.0,
            K_eta = 1.0,
            K_sigma = 1.0,
            K_kappa = 1.0,
            K_a0 = 1.0
        )
        
        densities = (
            u = 0.1,
            d = 0.1,
            s = 0.05,
            ubar = 0.1,
            dbar = 0.1,
            sbar = 0.05
        )
        
        # Test case 1: Low temperature
        q1 = QuarkParams((u=1.5, d=1.5, s=3.0), (u=0.2, d=0.2, s=0.2))
        t1 = ThermoParams(0.10, 0.3, 0.3, 0.0)
        
        result1 = Main.RelaxationTime.relaxation_times(
            q1, t1, K_coeffs;
            densities = densities
        )
        
        @test haskey(result1, :tau)
        @test all(isfinite, values(result1.tau))
        @test all(>(0), values(result1.tau))
        
        # Test case 2: High temperature
        q2 = QuarkParams((u=1.5, d=1.5, s=3.0), (u=0.4, d=0.4, s=0.4))
        t2 = ThermoParams(0.20, 0.7, 0.7, 0.0)
        
        result2 = Main.RelaxationTime.relaxation_times(
            q2, t2, K_coeffs;
            densities = densities
        )
        
        @test haskey(result2, :tau)
        @test all(isfinite, values(result2.tau))
        @test all(>(0), values(result2.tau))
        
        # Test case 3: Different quark masses
        q3 = QuarkParams((u=2.0, d=2.0, s=4.0), (u=0.3, d=0.3, s=0.3))
        t3 = ThermoParams(0.15, 0.5, 0.5, 0.0)
        
        result3 = Main.RelaxationTime.relaxation_times(
            q3, t3, K_coeffs;
            densities = densities
        )
        
        @test haskey(result3, :tau)
        @test all(isfinite, values(result3.tau))
        @test all(>(0), values(result3.tau))
    end
    
    @testset "Intermediate results consistency" begin
        # Test that intermediate functions in the chain work correctly with structs
        q_struct = QuarkParams(
            (u=1.52, d=1.52, s=3.04),  # m
            (u=0.3, d=0.3, s=0.3)       # μ
        )
        t_struct = ThermoParams(0.15, 0.5, 0.5, 0.0)
        
        K_coeffs = (
            K_pi = 1.0,
            K_K = 1.0,
            K_eta = 1.0,
            K_sigma = 1.0,
            K_kappa = 1.0,
            K_a0 = 1.0
        )
        
        # Ensure quark_params has A field for propagator calculations
        q_with_A = Main.RelaxationTime.ensure_quark_params_has_A(q_struct, t_struct)
        
        # Test AverageScatteringRate with struct parameters
        process = :us_to_us
        avg_rate = Main.RelaxationTime.AverageScatteringRate.average_scattering_rate(
            process, q_struct, t_struct, K_coeffs
        )
        
        @test isfinite(avg_rate)
        @test avg_rate >= 0
        
        # Test TotalCrossSection with struct parameters
        s = 5.0  # Mandelstam variable
        xs = Main.RelaxationTime.TotalCrossSection.total_cross_section(
            process, s, q_struct, t_struct, K_coeffs
        )
        
        @test isfinite(xs)
        @test xs >= 0
        
        # Test ScatteringAmplitude with struct parameters (needs A field)
        t_var = -1.0  # Mandelstam variable t
        amp_sq = Main.RelaxationTime.ScatteringAmplitude.scattering_amplitude_squared(
            process, s, t_var, q_with_A, t_struct, K_coeffs
        )
        
        @test isfinite(amp_sq)
        @test amp_sq >= 0
        
        # Test DifferentialCrossSection with struct parameters
        dxs = Main.RelaxationTime.DifferentialCrossSection.differential_cross_section(
            process, s, t_var, q_struct, t_struct, K_coeffs
        )
        
        @test isfinite(dxs)
        @test dxs >= 0
    end
    
    @testset "Mixed struct/NamedTuple usage" begin
        # Test mixing struct and NamedTuple parameters
        q_struct = QuarkParams(
            (u=1.52, d=1.52, s=3.04),  # m
            (u=0.3, d=0.3, s=0.3)       # μ
        )
        t_nt = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
        
        K_coeffs = (
            K_pi = 1.0,
            K_K = 1.0,
            K_eta = 1.0,
            K_sigma = 1.0,
            K_kappa = 1.0,
            K_a0 = 1.0
        )
        
        densities = (
            u = 0.1,
            d = 0.1,
            s = 0.05,
            ubar = 0.1,
            dbar = 0.1,
            sbar = 0.05
        )
        
        # Call with struct quark_params and NamedTuple thermo_params
        result1 = Main.RelaxationTime.relaxation_times(
            q_struct, t_nt, K_coeffs;
            densities = densities
        )
        
        @test haskey(result1, :tau)
        @test all(isfinite, values(result1.tau))
        @test all(>(0), values(result1.tau))
        
        # Call with NamedTuple quark_params and struct thermo_params
        q_nt = as_namedtuple(q_struct)
        t_struct = ThermoParams(0.15, 0.5, 0.5, 0.0)
        
        result2 = Main.RelaxationTime.relaxation_times(
            q_nt, t_struct, K_coeffs;
            densities = densities
        )
        
        @test haskey(result2, :tau)
        @test all(isfinite, values(result2.tau))
        @test all(>(0), values(result2.tau))
        
        # Verify both approaches give same results
        @test approx_equal(result1.tau.u, result2.tau.u)
        @test approx_equal(result1.tau.d, result2.tau.d)
        @test approx_equal(result1.tau.s, result2.tau.s)
    end
end
