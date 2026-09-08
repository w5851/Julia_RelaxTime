using Test
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","causal_gbu_infinite_profile.jl"))
const IP=CausalGBUInfiniteProfile

@testset "Exact monotone endpoints under floating affine cancellation" begin
    for a in (-12.,0.,1.4,3.832563530306745),b in (3.832563530306745,8.,36.),n in (1,512,1024)
        a<b || continue
        xs=IP.panel_nodes(a,b,n)
        @test first(xs)==a
        @test last(xs)==b
        @test all(diff(xs).>0)
    end
    a=3.832563530306745;b=nextfloat(a)
    @test IP.panel_nodes(a,b,1024)==[a,b]
end

@testset "Analytic square-root Cauchy window" begin
    for H in (0.3,1.,2.),z in (-3.,0.,0.15,0.7,3.,0.4im,0.6+0.3im,100+3im)
        r(x)=sqrt(x)*(1-x/H)
        direct=IP.I.C.mapped_cauchy(r,[0.,H],z;nodes=128)*pi
        @test IP.threshold_transform(z,H)≈direct atol=1e-11
    end
    @test IP.threshold_transform(1.,1.)≈-2/3 atol=1e-14
    for z in (0.4im,-3+0.5im,3+0.5im)
        t=(U=2.,H=0.5,positive=0.7,negative=-0.4)
        direct=IP.I.C.mapped_cauchy(x->IP.threshold_cut(t,x),[-2.5,-2.,2.,2.5],z;nodes=128)
        @test IP.threshold_value(t,z)≈direct atol=1e-11
    end
end

@testset "Infinite profile retains tails and vacuum jumps" begin
    bg=(m=(u=1.1,d=1.3,s=1.9),mu=(u=0.2,d=0.3,s=0.1),T=0.8,Phi=0.3,PhiBar=0.4,vacuum=2.5)
    for q in (0.,0.6)
        k=IP.I.kernel(bg,:K_plus,q)
        p=IP.profile(k;mesh=256)
        for x in (k.threshold+1e-5,k.threshold+1e-6)
            @test k.rho(x)/sqrt(x-k.threshold)≈p.threshold.positive rtol=1e-4
        end
        for z in (0.8im,2+0.8im,k.threshold, k.threshold+0.01, (k.threshold+k.landau)/2)
            @test IP.polarization(p,z)≈IP.I.polarization(k,z;nodes=128) atol=1e-5
        end
        @test any(!iszero,p.tailpositive)
        if q==0
            edge=hypot(bg.m.u,bg.vacuum)+hypot(bg.m.s,bg.vacuum)
            @test_throws ArgumentError IP.polarization(p,edge)
        end
    end
end
