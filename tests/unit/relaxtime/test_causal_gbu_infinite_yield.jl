using Test
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","causal_gbu_infinite_yield.jl"))
const testcausalgbuinfiniteyield_IY=CausalGBUInfiniteYield

@testset "Roots counted without a phase or unwrap" begin
    for sign in (-1.,1.)
        roots=testcausalgbuinfiniteyield_IY.bracket_roots(x->sign*(x-0.31)*(x-0.73),0.,1.)
        @test length(roots)==2
        @test roots≈[0.31,0.73] atol=1e-12
    end
    @test isempty(testcausalgbuinfiniteyield_IY.bracket_roots(x->x^2+1,0.,1.))
    @test_throws ArgumentError testcausalgbuinfiniteyield_IY.bracket_roots(identity,1.,0.)
    # Large finite Landau interval must not lose its thermal boundary layer.
    for h in (10.,1e3,1e8)
        @test testcausalgbuinfiniteyield_IY.I.thermal_range(x->exp(-x),0.,h,1.;nodes=96)≈1-exp(-h) atol=1e-12
    end
end

@testset "GBU bound/continuum integration-by-parts bookkeeping" begin
    # Unit step followed by a decreasing continuum, with one explicit root.
    T=0.8; mass=1.2; threshold=2.; endcut=4.
    delta(w)=pi*(endcut-w)/(endcut-threshold)
    g(w)=inv(expm1(w/T))
    direct=real(testcausalgbuinfiniteyield_IY.I.C.mapped_integral([threshold,endcut];nodes=96) do w
        g(w)*(1+g(w))/T*testcausalgbuinfiniteyield_IY.O.gbu_weight(delta(w))/pi
    end)-g(threshold)
    derivative=real(testcausalgbuinfiniteyield_IY.I.C.mapped_integral([threshold,endcut];nodes=96) do w
        -2g(w)*sin(delta(w))^2/(endcut-threshold)
    end)
    @test direct≈derivative atol=1e-12
    @test direct<0
    @test g(mass)+direct>0
end
