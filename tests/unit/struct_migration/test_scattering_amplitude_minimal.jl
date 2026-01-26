"""
Minimal test for ScatteringAmplitude struct equivalence - for quick verification.
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

# Load test utilities
include("test_utils.jl")
using .Main: QuarkParams, ThermoParams, as_namedtuple

@testset "ScatteringAmplitude Minimal Test" begin
    println("\n" * "="^70)
    println("Minimal Test: Scattering Amplitude Struct-NamedTuple Equivalence")
    println("="^70)
    
    # Fixed parameters for testing
    m_u = 1.0
    m_s = 3.0
    μ_u = 0.2
    μ_s = 0.2
    T = 0.15
    Φ = 0.5
    Φbar = 0.5
    s = 8.0
    t = -0.5
    
    # Compute A functions
    nodes_p, weights_p = gauleg(0.0, 20.0, 64)
    A_u = A(m_u, μ_u, T, Φ, Φbar, nodes_p, weights_p)
    A_s = A(m_s, μ_s, T, Φ, Φbar, nodes_p, weights_p)
    
    # Compute G functions and K_coeffs
    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)
    K_coeffs = calculate_effective_couplings(
        Constants_PNJL.G_fm2, Constants_PNJL.K_fm5, G_u, G_s
    )
    
    # Create parameters with A field
    q_nt_with_A = (m=(u=m_u, d=m_u, s=m_s), μ=(u=μ_u, d=μ_u, s=μ_s), A=(u=A_u, d=A_u, s=A_s))
    t_nt = (T=T, Φ=Φ, Φbar=Φbar, ξ=0.0)
    
    # Create struct parameters and add A field
    q_struct = QuarkParams((u=m_u, d=m_u, s=m_s), (u=μ_u, d=μ_u, s=μ_s))
    t_struct = ThermoParams(T, Φ, Φbar, 0.0)
    q_struct_with_A = merge(as_namedtuple(q_struct), (A=(u=A_u, d=A_u, s=A_s),))
    
    # Test uu_to_uu process
    process = :uu_to_uu
    
    M_squared_struct = scattering_amplitude_squared(
        process, s, t, q_struct_with_A, t_nt, K_coeffs
    )
    M_squared_nt = scattering_amplitude_squared(
        process, s, t, q_nt_with_A, t_nt, K_coeffs
    )
    
    println("Process: $process")
    println("  Struct result:     $M_squared_struct")
    println("  NamedTuple result: $M_squared_nt")
    println("  Difference:        $(abs(M_squared_struct - M_squared_nt))")
    
    @test isapprox(M_squared_struct, M_squared_nt, rtol=1e-12, atol=1e-14)
    @test M_squared_struct >= 0.0
    @test isfinite(M_squared_struct)
    
    # Test uubar_to_uubar process
    process = :uubar_to_uubar
    
    M_squared_struct = scattering_amplitude_squared(
        process, s, t, q_struct_with_A, t_nt, K_coeffs
    )
    M_squared_nt = scattering_amplitude_squared(
        process, s, t, q_nt_with_A, t_nt, K_coeffs
    )
    
    println("\nProcess: $process")
    println("  Struct result:     $M_squared_struct")
    println("  NamedTuple result: $M_squared_nt")
    println("  Difference:        $(abs(M_squared_struct - M_squared_nt))")
    
    @test isapprox(M_squared_struct, M_squared_nt, rtol=1e-12, atol=1e-14)
    @test M_squared_struct >= 0.0
    @test isfinite(M_squared_struct)
    
    println("\n✓ All minimal tests passed!")
    println("="^70)
end
