using Test
const _PEREIRA_ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
include(joinpath(_PEREIRA_ROOT,"scripts","analysis","relaxtime","causal_gbu_pereira_reference.jl"))
const testcausalgbupereirareference_PR=CausalGBUPereiraReference
for (name,path) in ((:Constants_PNJL,"constants/Constants_PNJL.jl"),
                   (:PNJLQuarkDistributions,"models/pnjl_physics/QuarkDistribution.jl"),
                   (:OneLoopIntegrals,"relaxtime/OneLoopIntegrals.jl"),
                   (:CausalSpectralBubble,"relaxtime/CausalSpectralBubble.jl"))
    isdefined(Main,name) || Base.include(Main,joinpath(_PEREIRA_ROOT,"src",path))
end
const testcausalgbupereirareference_PCB=Main.CausalSpectralBubble

@testset "Pereira: occupation reduction and exact pair/scattering algebra" begin
    for x in (-1000.,-10.,-0.2,0.,0.2,10.,1000.)
        fd=x>=0 ? exp(-x)/(1+exp(-x)) : 1/(1+exp(x))
        @test testcausalgbupereirareference_PR.occupation(x,1.,1.)≈fd atol=2e-16
        @test testcausalgbupereirareference_PR.occupation(x,0.,0.)≈(x>=0 ? exp(-3x)/(1+exp(-3x)) : 1/(1+exp(3x))) atol=2e-16
    end
    for (phi,bar) in ((1.,1.),(0.3,0.4)),mu in (-0.3,0.2)
        n=testcausalgbupereirareference_PR.occupations(1.1,mu,0.8,phi,bar)
        @test n.quark≈Main.PNJLQuarkDistributions.quark_distribution(1.1,mu,0.8,phi,bar) atol=1e-15
        @test n.anti≈Main.PNJLQuarkDistributions.antiquark_distribution(1.1,mu,0.8,phi,bar) atol=1e-15
    end
    E1,E2,z=1.3,1.8,0.7+0.4im
    for c in (0.,1.)
        n1=(quark=0.13,anti=0.09); n2=(quark=0.07,anti=0.17)
        ref=sum(testcausalgbupereirareference_PR.pair_scattering(E1,E2,z,n1,n2,c))
        n1s=(c-n1.anti,n1.quark); n2s=(c-n2.anti,n2.quark)
        residue=sum(s*t*(n2s[j]-n1s[i])/(z-s*E1+t*E2)
                    for (i,s) in enumerate((-1,1)),(j,t) in enumerate((-1,1)))
        @test ref≈residue atol=1e-15
        W1,W2=c-n1.anti-n2.quark,c-n1.quark-n2.anti
        Dq,Da=n1.quark-n2.quark,n1.anti-n2.anti
        @test W1-W2-Dq+Da≈0 atol=1e-15 # no 1/z term in B0
        moment=-(E1+E2)*(W1+W2)+(E2-E1)*(Dq+Da)
        contact=-2*E2*(c-n1.quark-n1.anti)-2*E1*(c-n2.quark-n2.anti)
        @test moment≈contact atol=2e-15 # fixes the constant cancellation in Pi
    end
    @test testcausalgbupereirareference_PR.quadratic_intervals(1.,0.,-1.,-2.,2.)==[(-1.,1.)]
    @test testcausalgbupereirareference_PR.quadratic_intervals(-1.,0.,1.,-2.,2.)==[(-2.,-1.),(1.,2.)]
    @test testcausalgbupereirareference_PR.quadratic_intervals(0.,1.,0.,-2.,2.)==[(-2.,-0.)]
    @test isempty(testcausalgbupereirareference_PR.quadratic_intervals(1.,0.,1.,-2.,2.))
end

@testset "Pereira: exact lens volume and same-domain single-line terms" begin
    for q in (0.,0.01,0.8,3.,5.9,6.,6.1)
        L=3.
        r=testcausalgbupereirareference_PR.cylindrical_reference(q,1.,0.2,1.4,-0.1,0.8,L,[0.4+0.6im];
                                  component=:vacuum,nz=64,ny=64)
        expected=q<2L ? pi*(2L-q)^2*(4L+q)/12 : 0.
        @test r.volume≈expected atol=1e-11
        @test r.A1≈testcausalgbupereirareference_PR.vacuum_A(q,1.,L) atol=2e-11
        @test r.A2≈testcausalgbupereirareference_PR.vacuum_A(q,1.4,L) atol=2e-11
        @test !r.production_authorized
    end
    r=testcausalgbupereirareference_PR.cylindrical_reference(0.8,1.,0.2,2.7,0.5,0.8,3.,[0.4+0.6im];component=:vacuum)
    @test r.A1≈testcausalgbupereirareference_PR.vacuum_A(0.8,1.,3.) atol=2e-11
    # The inherited single-line term is not a q-independent equilibrium A.
    @test abs(r.A1-testcausalgbupereirareference_PR.vacuum_A(0.,1.,3.))>1.
    @test_throws ArgumentError testcausalgbupereirareference_PR.cylindrical_reference(-1.,1.,0.,1.,0.,0.8,3.,[1im])
    @test_throws ArgumentError testcausalgbupereirareference_PR.cylindrical_reference(1.,1.,0.,1.,0.,0.8,3.,[1.])
    @test_throws ArgumentError testcausalgbupereirareference_PR.vacuum_A(1.,0.,3.)
end

@testset "Pereira: independent complex loop, P/S contact, and PNJL extension" begin
    zs=[-1.1+0.7im,0.4+0.6im,3.2+0.8im]
    for q in (0.,0.8,3.),(phi,bar) in ((1.,1.),(0.3,0.4))
        a,u,b,v,T,L=1.,0.2,1.4,-0.1,0.8,3.
        r=testcausalgbupereirareference_PR.cylindrical_reference(q,a,u,b,v,T,L,zs;Phi=phi,PhiBar=bar,nz=80,ny=80)
        rev=testcausalgbupereirareference_PR.cylindrical_reference(q,b,v,a,u,T,L,-conj.(zs);Phi=phi,PhiBar=bar,nz=80,ny=80)
        @test r.b0≈conj.(rev.b0) atol=2e-10
        for channel in (:P,:S)
            g=testcausalgbupereirareference_PCB.build_spectral_bubble(q,a,u,b,v,T;Phi=phi,PhiBar=bar,channel=channel,
                vacuum_cutoff_inv_fm=L,thermal_cutoff_inv_fm=L,momentum_nodes=96,angle_nodes=64)
            @test g.contact_inv_fm2≈r.contact atol=2e-10
            for (i,z) in enumerate(zs)
                d=testcausalgbupereirareference_PCB.spectral_bubble(g,real(z)-u+v;eta_inv_fm=imag(z))
                @test d.B0≈r.b0[i] atol=2e-9
                @test d.value≈(channel===:P ? r.pi_p[i] : r.pi_s[i]) atol=2e-9
            end
        end
    end
end

@testset "Pereira: C18/C19 retarded cuts, support, and thermal decomposition" begin
    for q in (0.,0.8,3.,6.),l in (-4.0,-0.2,0.,0.2,4.0),component in (:full,:vacuum,:thermal)
        a,u,b,v,T,L,phi,bar=1.,0.2,1.4,-0.1,0.8,3.,0.3,0.4
        r=testcausalgbupereirareference_PR.reference_cut(l,q,a,u,b,v,T,L;Phi=phi,PhiBar=bar,component=component)
        c=Main.OneLoopIntegrals.B0_spectral_cut(l,q,a,u,b,v,T;Φ=phi,Φbar=bar,
            pmax_inv_fm=L,component=component)
        @test r.pair≈c.pair atol=2e-11
        @test r.landau≈c.landau atol=2e-11
        @test !r.production_authorized
    end
    for q in (0.,0.8),l in (0.2,3.0)
        f(c)=testcausalgbupereirareference_PR.reference_cut(l,q,1.,0.2,1.4,-0.1,0.8,3.;component=c).imaginary
        @test f(:full)≈f(:vacuum)+f(:thermal) atol=1e-13
        plus=testcausalgbupereirareference_PR.reference_cut(3.,q,1.,0.,1.4,0.,0.8,3.;component=:vacuum)
        minus=testcausalgbupereirareference_PR.reference_cut(-3.,q,1.,0.,1.4,0.,0.8,3.;component=:vacuum)
        @test plus.pair>0 && minus.pair<0 # mass-i0 would instead be even
        @test plus.pair≈-minus.pair atol=1e-13
    end
    @test testcausalgbupereirareference_PR.reference_cut(0.,0.,1.,0.,1.,0.,0.8,3.).static_degeneracy_unresolved
    @test_throws ArgumentError testcausalgbupereirareference_PR.reference_cut(NaN,1.,1.,0.,1.,0.,0.8,3.)
end

# Including the executable also catches parse/load regressions without running main.
include(joinpath(_PEREIRA_ROOT,"scripts","analysis","relaxtime","audit_causal_gbu_pereira.jl"))
const testcausalgbupereirareference_PA=CausalGBUPereiraAudit
@testset "Pereira: runner validation and synthetic adapters" begin
    @test_throws ArgumentError testcausalgbupereirareference_PA.main(nodes=15)
    @test_throws ArgumentError testcausalgbupereirareference_PA.main(qs=(NaN,))
    @test_throws ArgumentError testcausalgbupereirareference_PA.main(qs=())
    @test_throws ErrorException testcausalgbupereirareference_PA.R.start_output(_PEREIRA_ROOT)
    bg=(m=(u=1.,d=1.4,s=1.8),mu=(u=0.2,d=-0.1,s=0.),
        T=0.8,vacuum=3.)
    zs=[0.4+0.6im,3.2+0.8im]
    combined=testcausalgbupereirareference_PA.reference(bg,:pi_plus,0.8,zs,0.3,0.4,3.,48)
    full=testcausalgbupereirareference_PR.cylindrical_reference(0.8,1.,0.2,1.4,-0.1,0.8,3.,zs;
        Phi=0.3,PhiBar=0.4,nz=48,ny=48)
    @test combined.b0≈full.b0 atol=2e-12
    @test combined.pi≈full.pi_p atol=2e-12
    @test combined.contact≈full.contact atol=2e-12
    for l in (-3.,-0.2,0.2,3.)
        c=testcausalgbupereirareference_PA.cut(bg,:pi_plus,0.8,l,0.3,0.4,3.,64)
        f=testcausalgbupereirareference_PR.reference_cut(l,0.8,1.,0.2,1.4,-0.1,0.8,3.;Phi=0.3,PhiBar=0.4)
        @test c.pair≈f.pair atol=2e-13
        @test c.landau≈f.landau atol=2e-13
    end
end
