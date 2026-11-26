using Test
using StaticArrays
using ForwardDiff
using SpecialFunctions: loggamma

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "magnetic thermal and small-mass limits" begin
    mi = Models.MagneticThermodynamics.MagneticIntegrals
    hbarc = 197.3269804
    eB = .4*(1000/hbarc)^2
    q = 2/3
    qB = q*eB

    # Analytic massless lowest-Landau-level Fermi gas (Phi=Phibar=1).
    # Density is exactly 3*qB*mu/(2pi^2), at any positive T.
    mu = 2.3
    for T in (.015,.5), mass in (1e-6,.001)
        rho = mi.density_flavor_landau(mass,mu,T,1.,1.,q,eB;
            n_max=0,p_num=64,pz_max=40.)
        expected = 3qB*mu/(2pi^2)
        @test isapprox(rho,expected;rtol=2e-7,atol=2e-9)
    end
    # Direct density and derivative of the integrated pressure, including
    # moving Fermi-surface quadrature nodes, must be the same observable.
    for mu in (-2.3,0.,2.3)
        f = z->mi.omegat_flavor_landau(.035,z,.025,.15,.23,q,eB;
            n_max=4,p_num=64,pz_max=40.)
        rho = mi.density_flavor_landau(.035,mu,.025,.15,.23,q,eB;
            n_max=4,p_num=64,pz_max=40.)
        @test isapprox(-ForwardDiff.derivative(f,mu),rho;rtol=1e-8,atol=2e-10)
    end
    # Independent log-gamma formula for the MFIR magnetic condensate.
    for x in (1e-8,7.60072e-5,.01,.9999,1.0001,10.)
        mass = sqrt(2qB*x)
        expected = -3mass*qB/(2pi^2)*(loggamma(x)-log(2pi)/2+x-(x-.5)*log(x))
        got = ForwardDiff.derivative(m->mi.omega_magnetic_mfir(m,q,eB;zeta_num=256),mass)
        @test isapprox(got,expected;rtol=2e-8,atol=2e-9)
    end
    @test_throws ArgumentError mi.omegat_flavor_landau(1.,1.,0.,.1,.1,q,eB;n_max=1)
    @test_throws ArgumentError mi.density_flavor_landau(1.,1.,.1,.1,.1,q,eB;n_max=-1)
end
Models.pnjl_module()
const PNJL = Models.pnjl_module()
const alpha_n = getproperty(PNJL, :alpha_n)
const energy_landau = getproperty(PNJL, :energy_landau)
const resolve_nmax_from_cutoff = getproperty(PNJL, :resolve_nmax_from_cutoff)
const smooth_cutoff = getproperty(PNJL, :smooth_cutoff)

@testset "magnetic energy levels" begin
    @test alpha_n(0) == 1.0
    @test alpha_n(1) == 2.0
    @test alpha_n(5) == 2.0

    m = 1.2
    pz = 0.8
    q = 2 / 3
    eB = 0.15

    E0 = energy_landau(m, pz, 0, q, eB)
    E1 = energy_landau(m, pz, 1, q, eB)
    @test isfinite(E0)
    @test isfinite(E1)
    @test E1 > E0

    nmax = resolve_nmax_from_cutoff(m, 1.0, q, eB)
    @test nmax >= 0

    val = smooth_cutoff(1.0; Λ=3.0, N=10)
    @test 0.0 < val <= 1.0
end
