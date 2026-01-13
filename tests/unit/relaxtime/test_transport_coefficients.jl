using Test

include("../../../src/relaxtime/TransportCoefficients.jl")
using .TransportCoefficients

const QUARK_PARAMS = (m=(u=0.3,d=0.3,s=0.5), μ=(u=0.2,d=0.2,s=0.2))
const THERMO_PARAMS = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)

const TAU_ZERO = (u=0.0,d=0.0,s=0.0,ubar=0.0,dbar=0.0,sbar=0.0)
const TAU_ONE = (u=1.0,d=1.0,s=1.0,ubar=1.0,dbar=1.0,sbar=1.0)

@testset "TransportCoefficients: eta/sigma basic" begin
    eta0 = shear_viscosity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ZERO, p_nodes=16, p_max=10.0)
    sigma0 = electric_conductivity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ZERO, p_nodes=16, p_max=10.0)
    @test eta0 == 0.0
    @test sigma0 == 0.0

    eta1 = shear_viscosity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, p_nodes=16, p_max=10.0)
    sigma1 = electric_conductivity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, p_nodes=16, p_max=10.0)
    @test isfinite(eta1)
    @test isfinite(sigma1)
    @test eta1 > 0.0
    @test sigma1 > 0.0
end

@testset "TransportCoefficients: sigma scales with q^2" begin
    charges1 = (u=2/3, d=-1/3, s=-1/3)
    charges2 = (u=4/3, d=-2/3, s=-2/3)

    σ1 = electric_conductivity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, charges=charges1, p_nodes=16, p_max=10.0)
    σ2 = electric_conductivity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, charges=charges2, p_nodes=16, p_max=10.0)

    # charges2 = 2 * charges1 => q^2 scales by 4
    @test isapprox(σ2 / σ1, 4.0; rtol=1e-3)
end

