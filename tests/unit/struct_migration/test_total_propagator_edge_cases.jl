"""
Unit tests for TotalPropagator module edge cases with struct parameters.

Tests:
- Different scattering processes and channels
- Different meson combinations
- Kinematic edge cases
- Error handling
"""

using Test

# Load required modules
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src/relaxtime"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../../src/utils"))

# Load ParticleSymbols first (required by TotalPropagator)
if !isdefined(Main, :ParticleSymbols)
    Base.include(Main, joinpath(@__DIR__, "../../../src/utils/ParticleSymbols.jl"))
end

include("../../../src/Constants_PNJL.jl")
include("../../../src/relaxtime/TotalPropagator.jl")
include("../../../src/relaxtime/EffectiveCouplings.jl")
include("../../../src/relaxtime/OneLoopIntegrals.jl")
include("../../../src/integration/GaussLegendre.jl")

using .TotalPropagator
using .Constants_PNJL

# Load test utilities
include("test_utils.jl")
using .Main: QuarkParams, ThermoParams, as_namedtuple, approx_equal

@testset "TotalPropagator Edge Cases" begin
    # Test fixtures
    q_nt = get_test_quark_params()
    t_nt = get_test_thermo_params()
    K_nt = get_test_K_coeffs()
    
    q_struct = QuarkParams(q_nt)
    t_struct = ThermoParams(t_nt)
    
    @testset "Different scattering processes" begin
        # Test various scattering processes
        processes = [:uu_to_uu, :ud_to_ud, :us_to_us, :dd_to_dd, :ds_to_ds, :ss_to_ss]
        channels = [:s, :t, :u]
        
        for process in processes
            for channel in channels
                # Test with struct parameters
                result_struct = calculate_all_propagators(
                    process, 1.0, 1.0,
                    q_struct, t_struct, K_nt
                )
                
                # Test with NamedTuple parameters
                result_nt = calculate_all_propagators(
                    process, 1.0, 1.0,
                    q_nt, t_nt, K_nt
                )
                
                # Verify both produce valid results
                @test result_struct isa NamedTuple
                @test result_nt isa NamedTuple
                
                # Verify equivalence
                @test approx_equal(result_struct, result_nt)
            end
        end
    end
    
    @testset "Different meson combinations" begin
        # Test different meson lists
        meson_lists = [
            [:pi],
            [:K],
            [:eta],
            [:pi, :K],
            [:pi, :eta],
            [:K, :eta],
            [:pi, :K, :eta]
        ]
        
        for mesons in meson_lists
            # Test with struct parameters
            result_struct = total_propagator_simple(
                :uu_to_uu, :s, mesons,
                1.0, 1.0,
                q_struct, t_struct, K_nt
            )
            
            # Test with NamedTuple parameters
            result_nt = total_propagator_simple(
                :uu_to_uu, :s, mesons,
                1.0, 1.0,
                q_nt, t_nt, K_nt
            )
            
            # Verify both produce valid results
            @test result_struct isa Float64
            @test result_nt isa Float64
            @test isfinite(result_struct)
            @test isfinite(result_nt)
            
            # Verify equivalence
            @test approx_equal(result_struct, result_nt)
        end
    end
    
    @testset "Kinematic edge cases" begin
        # Test near-zero momentum
        @testset "Near-zero momentum" begin
            k0_small = 0.01
            k_small = 0.01
            
            result_struct = calculate_all_propagators(
                :uu_to_uu, k0_small, k_small,
                q_struct, t_struct, K_nt
            )
            
            result_nt = calculate_all_propagators(
                :uu_to_uu, k0_small, k_small,
                q_nt, t_nt, K_nt
            )
            
            @test result_struct isa NamedTuple
            @test result_nt isa NamedTuple
            @test approx_equal(result_struct, result_nt)
        end
        
        # Test large momentum
        @testset "Large momentum" begin
            k0_large = 10.0
            k_large = 10.0
            
            result_struct = calculate_all_propagators(
                :uu_to_uu, k0_large, k_large,
                q_struct, t_struct, K_nt
            )
            
            result_nt = calculate_all_propagators(
                :uu_to_uu, k0_large, k_large,
                q_nt, t_nt, K_nt
            )
            
            @test result_struct isa NamedTuple
            @test result_nt isa NamedTuple
            @test approx_equal(result_struct, result_nt)
        end
        
        # Test different k0/k ratios
        @testset "Different k0/k ratios" begin
            test_cases = [
                (k0=1.0, k=0.5),   # k0 > k
                (k0=0.5, k=1.0),   # k0 < k
                (k0=1.0, k=1.0),   # k0 = k
            ]
            
            for (k0, k) in test_cases
                result_struct = calculate_all_propagators(
                    :uu_to_uu, k0, k,
                    q_struct, t_struct, K_nt
                )
                
                result_nt = calculate_all_propagators(
                    :uu_to_uu, k0, k,
                    q_nt, t_nt, K_nt
                )
                
                @test approx_equal(result_struct, result_nt)
            end
        end
    end
    
    @testset "CMS momentum calculation" begin
        # Test CMS momentum for different processes and channels
        processes = [:uu_to_uu, :ud_to_ud, :us_to_us]
        channels = [:s, :t, :u]
        
        s_test = 31.0
        t_test = -5.0
        
        for process in processes
            for channel in channels
                # Test with struct parameters
                result_struct = calculate_cms_momentum(
                    process, s_test, t_test, channel,
                    q_struct
                )
                
                # Test with NamedTuple parameters
                result_nt = calculate_cms_momentum(
                    process, s_test, t_test, channel,
                    q_nt
                )
                
                # Verify both produce valid results
                @test result_struct isa NamedTuple
                @test result_nt isa NamedTuple
                @test haskey(result_struct, :k0)
                @test haskey(result_struct, :k_norm)
                
                # Verify equivalence
                @test approx_equal(result_struct.k0, result_nt.k0)
                @test approx_equal(result_struct.k_norm, result_nt.k_norm)
            end
        end
    end
    
    @testset "Mixed propagator" begin
        # Test mixed propagator with different meson channels
        meson_channels = [:pi, :K, :eta]
        
        for meson_channel in meson_channels
            result_struct = total_propagator_mixed(
                :uu_to_uu, :s, meson_channel,
                1.0, 1.0,
                q_struct, t_struct, K_nt
            )
            
            result_nt = total_propagator_mixed(
                :uu_to_uu, :s, meson_channel,
                1.0, 1.0,
                q_nt, t_nt, K_nt
            )
            
            @test result_struct isa Float64
            @test result_nt isa Float64
            @test isfinite(result_struct)
            @test isfinite(result_nt)
            @test approx_equal(result_struct, result_nt)
        end
    end
    
    @testset "Auto propagator" begin
        # Test auto propagator (automatically selects mesons)
        processes = [:uu_to_uu, :ud_to_ud, :us_to_us]
        channels = [:s, :t, :u]
        
        for process in processes
            for channel in channels
                result_struct = total_propagator_auto(
                    process, channel,
                    1.0, 1.0,
                    q_struct, t_struct, K_nt
                )
                
                result_nt = total_propagator_auto(
                    process, channel,
                    1.0, 1.0,
                    q_nt, t_nt, K_nt
                )
                
                @test result_struct isa Float64
                @test result_nt isa Float64
                @test isfinite(result_struct)
                @test isfinite(result_nt)
                @test approx_equal(result_struct, result_nt)
            end
        end
    end
    
    @testset "Flavor factor" begin
        # Test flavor factor extraction
        processes = [:uu_to_uu, :ud_to_ud, :us_to_us, :dd_to_dd, :ds_to_ds, :ss_to_ss]
        
        for process in processes
            # Flavor factor should be the same for struct and NamedTuple
            factor = get_flavor_factor(process)
            
            @test factor isa Float64
            @test factor > 0.0
            @test isfinite(factor)
        end
    end
    
    @testset "All propagators by channel" begin
        # Test calculating all propagators organized by channel
        result_struct = calculate_all_propagators_by_channel(
            :uu_to_uu, 1.0, 1.0,
            q_struct, t_struct, K_nt
        )
        
        result_nt = calculate_all_propagators_by_channel(
            :uu_to_uu, 1.0, 1.0,
            q_nt, t_nt, K_nt
        )
        
        @test result_struct isa NamedTuple
        @test result_nt isa NamedTuple
        
        # Should have s, t, u channels
        @test haskey(result_struct, :s)
        @test haskey(result_struct, :t)
        @test haskey(result_struct, :u)
        
        # Verify equivalence
        @test approx_equal(result_struct, result_nt)
    end
end
