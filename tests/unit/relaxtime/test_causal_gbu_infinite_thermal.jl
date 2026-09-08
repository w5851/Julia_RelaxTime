using Test
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","causal_gbu_infinite_thermal.jl"))
const testcausalgbuinfinitethermal_IT=CausalGBUInfiniteThermal

@testset "Infinite-domain quadrature and exact domains" begin
    @test testcausalgbuinfinitethermal_IT.integrate_range(x->exp(-x),0.,Inf;nodes=96)≈1 atol=1e-12
    @test testcausalgbuinfinitethermal_IT.integrate_range(x->x^2*exp(-x),2.,Inf;nodes=96)≈10exp(-2) atol=1e-12
    @test testcausalgbuinfinitethermal_IT.integrate_range(x->x^2,0.,2.;nodes=32)≈8/3 atol=1e-14
    @test testcausalgbuinfinitethermal_IT.infinite_intervals(-1.,0.,1.,0.)==[(1.,Inf)]
    @test testcausalgbuinfinitethermal_IT.infinite_intervals(1.,-3.,2.,0.)==[(1.,2.)]
    @test testcausalgbuinfinitethermal_IT.infinite_intervals(0.,-1.,2.,0.)==[(2.,Inf)]
    @test testcausalgbuinfinitethermal_IT.infinite_intervals(0.,1.,-2.,0.)==[(0.,2.)]
    @test isempty(testcausalgbuinfinitethermal_IT.infinite_intervals(1.,0.,1.,0.))
    @test_throws ArgumentError testcausalgbuinfinitethermal_IT.integrate_range(identity,0.,1.;scale=0.)
end

@testset "Independent infinite radial response and signed-zero panels" begin
    bg=(m=(u=1.1,d=1.3,s=1.9),mu=(u=0.2,d=0.3,s=0.1),T=0.8,Phi=0.3,PhiBar=0.4,vacuum=2.5)
    for q in (0.,0.6),ch in (:pi_plus,:K_plus)
        k=testcausalgbuinfinitethermal_IT.kernel(bg,ch,q;cut_nodes=192)
        @test all(diff(k.edges).>0)
        for z in (0.8im,2+0.8im)
            p=testcausalgbuinfinitethermal_IT.polarization(k,z;nodes=64)
            radial=testcausalgbuinfinitethermal_IT.radial_polarization(bg,ch,q,z;nodes=192)
            @test p≈radial atol=1e-8
            @test testcausalgbuinfinitethermal_IT.polarization(k,z;nodes=96)≈p atol=1e-8
        end
    end
    @test_throws ArgumentError testcausalgbuinfinitethermal_IT.radial_polarization(bg,:pi_plus,0.,1.)
end

@testset "Infinite thermal cut equals large-L reference away from its artificial edge" begin
    bg=(m=(u=1.1,d=1.3,s=1.9),mu=(u=0.2,d=0.3,s=0.1),T=0.8,Phi=0.3,PhiBar=0.4,vacuum=2.5)
    for ch in testcausalgbuinfinitethermal_IT.R.CHANNELS,q in (0.,0.6),x in (-5.,-0.4,0.2,3.5,7.)
        actual=testcausalgbuinfinitethermal_IT.total_cut(bg,ch,q,x;nodes=96)
        reference=testcausalgbuinfinitethermal_IT.A.direct_cut(bg,ch,q,x,32.;nodes=128).imaginary
        @test actual≈reference atol=1e-10
        @test testcausalgbuinfinitethermal_IT.total_cut(bg,ch,q,BigFloat(x);nodes=96)≈actual atol=1e-12
    end
    # No edge at a nominal Lth=10: the exact thermal tail is nonzero beyond it.
    @test testcausalgbuinfinitethermal_IT.thermal_cut(bg,:pi_plus,0.,22.)<0
    for (ch,anti) in ((:pi_plus,:pi_minus),(:K_plus,:K_minus)),q in (0.,0.6),x in (-5.,-0.4,0.2,3.5,7.)
        @test testcausalgbuinfinitethermal_IT.total_cut(bg,ch,q,x;nodes=96)≈-testcausalgbuinfinitethermal_IT.total_cut(bg,anti,q,-x;nodes=96) atol=1e-12
    end
    b=testcausalgbuinfinitethermal_IT.phase_band_bound(20.,21.,0.8,1.;branch_bound_verified=true)
    @test 0<b.shell_absolute_bound_inv_fm2<1e-12
    @test b.bound_verified && !b.production_authorized
    @test !testcausalgbuinfinitethermal_IT.phase_band_bound(20.,21.,0.8,1.).bound_verified
end
