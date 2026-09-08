using Test
const _RC_ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
include(joinpath(_RC_ROOT,"scripts","analysis","relaxtime","causal_gbu_regulator_checks.jl"))
const testcausalgburegulatorchecks_RC=CausalGBURegulatorChecks

@testset "Regulator geometry and exact cut restriction" begin
    for (a,b) in ((1.0,1.4),(1.4,1.0),(1.0,1.0)),q in (0.0,0.5,3.0,7.0)
        L=3.0
        sp=testcausalgburegulatorchecks_RC.domain_support(q,a,b,L,:centered)
        for p in range(0,L;length=15),x in range(-1,1;length=15)
            e1=sqrt(p^2+q^2/4+p*q*x+a^2)
            e2=sqrt(p^2+q^2/4-p*q*x+b^2)
            @test sp.pair[1]-1e-12<=e1+e2<=sp.pair[2]+1e-12
            @test sp.landau[1]-1e-12<=e1-e2<=sp.landau[2]+1e-12
        end
        for l in (0.2,1.5,4.0,8.0),s in (-1,1),t in (-1,1)
            q==0 && continue
            for (u,v) in testcausalgburegulatorchecks_RC.centered_intervals(l,q,a,b,L,s,t)
                for E1 in range(u+1e-10,v-1e-10;length=5)
                    E2=t*(s*E1-l)
                    @test E1>=a && E2>=b
                    p2=E1^2-a^2
                    r2=E2^2-b^2
                    P2=(p2+r2)/2-q^2/4
                    @test -1e-10<=P2<=L^2+1e-10
                    @test abs((p2+q^2-r2)/(2sqrt(p2)*q))<=1+1e-9
                end
            end
        end
    end
    @test testcausalgburegulatorchecks_RC.domain_support(6.0,1.0,1.4,3.0,:two_line)===nothing
    @test testcausalgburegulatorchecks_RC.merge_intervals([(1.,2.),(-2.,-1.),(1.5,3.)])==[(-2.,-1.),(1.,3.)]
    @test testcausalgburegulatorchecks_RC.analytic_gaps([(-3.,-2.),(2.,3.)],0.5;upper=5.0)==[(0.,1.5),(2.5,5.)]
    @test_throws ArgumentError testcausalgburegulatorchecks_RC.domain_support(1.,1.,1.4,3.,:invalid)
    @test_throws ArgumentError testcausalgburegulatorchecks_RC.domain_support(NaN,1.,1.4,3.,:centered)
    @test_throws ArgumentError testcausalgburegulatorchecks_RC.centered_intervals(NaN,1.,1.,1.4,3.,1,-1)
    @test_throws ArgumentError testcausalgburegulatorchecks_RC.centered_cut(0.3,NaN,1.,0.2,1.4,-0.1,0.8,0.3,0.4,3.)
    @test_throws ArgumentError testcausalgburegulatorchecks_RC.centered_atoms(0.8,1.,0.2,1.4,-0.1,0.8,0.3,0.4,3.,NaN)
end

@testset "Full support gaps are counted independently of a normal gap" begin
    s=testcausalgburegulatorchecks_RC.R.Settings()
    cuts=[(-3.,-2.),(2.,3.)]
    bubble=(inverse=z->(z-0.5)*(z-2.0)*(z-4.0),grid=(q_inv_fm=0.,),shift=0.5)
    rows=testcausalgburegulatorchecks_RC.count_support_gaps(bubble,s,cuts;upper=8.0)
    @test length(rows)==2
    @test all(r.passed for r in rows)
    @test sum(r.root_count for r in rows)==2 # the zero in the cut is NOT a bound-state claim
    even=(inverse=z->(z-0.54321)^2,grid=(q_inv_fm=0.,),shift=0.5)
    miss=testcausalgburegulatorchecks_RC.count_support_gaps(even,s,cuts;upper=8.0)
    @test !first(miss).passed # simple real-root scanner cannot certify a double zero
    @test first(miss).contour_count==2
    near=(inverse=z->z-1.5e-6,grid=(q_inv_fm=0.,),shift=0.5)
    boundary=testcausalgburegulatorchecks_RC.count_support_gaps(near,s,cuts;upper=8.0)
    @test first(boundary).root_count==1 # margin must not be applied twice
    @test first(boundary).passed
end

@testset "Centered cut Cauchy reconstruction versus independent loop" begin
    bg=(m=(u=1.0,d=1.1,s=1.4),mu=(u=0.15,d=0.1,s=0.05),T=0.8,Phi=0.3,PhiBar=0.4,
        coupling=Dict(c=>0.20 for c in testcausalgburegulatorchecks_RC.R.CHANNELS),vacuum=3.0)
    for q in (0.0,0.8)
        # q=0 has a hard-cutoff jump; mesh64 gave 5.3--6.0e-4 error.
        # Resolve the endpoint cell, preserving the existing error target.
        low=testcausalgburegulatorchecks_RC.centered_profile(bg,:K_plus,q;mesh=64,ne=64)
        high=testcausalgburegulatorchecks_RC.centered_profile(bg,:K_plus,q;mesh=128,ne=64)
        g=testcausalgburegulatorchecks_RC.centered_atoms(q,1.,0.15,1.4,0.05,0.8,0.3,0.4,3.,24.;np=192,nx=96)
        for z in (0.5+0.6im,3.0+0.6im)
            target=testcausalgburegulatorchecks_RC.atom_value(g,z).value
            coarse=abs(testcausalgburegulatorchecks_RC.R.cauchy_transform(low.profile,z)-target)
            fine=abs(testcausalgburegulatorchecks_RC.R.cauchy_transform(high.profile,z)-target)
            @test fine<coarse
            @test fine<5e-4
        end
        for (left,right) in testcausalgburegulatorchecks_RC.analytic_gaps(high.support,high.shift)
            for w in range(left+1e-5,right-1e-5;length=5)
                @test imag(high.inverse(w))==0
            end
        end
    end
    # Exact Poisson smoothing identity, checked by independent real integration.
    p=testcausalgburegulatorchecks_RC.R.PiecewiseSpectralFunction([-1.,0.,1.],[0.,1.,0.])
    for gamma in (0.2,0.7,1.0),center in (0.0,0.4)
        integral=0.0
        for (left,right) in ((-1.,0.),(0.,1.))
            xs,ws=testcausalgburegulatorchecks_RC.R.gauleg(left,right,128)
            integral+=sum(ws[i]*(1-abs(xs[i]))*gamma/(pi*((xs[i]-center)^2+gamma^2)) for i in eachindex(xs))
        end
        @test integral≈imag(testcausalgburegulatorchecks_RC.R.cauchy_transform(p,center+gamma*im)) atol=1e-13
    end
end

@testset "Independent centered causal loop and contact identity" begin
    for q in (0.0,0.8)
        g=testcausalgburegulatorchecks_RC.centered_atoms(q,1.,0.2,1.4,-0.1,0.8,0.3,0.4,3.,6.;np=64,nx=48)
        rev=testcausalgburegulatorchecks_RC.centered_atoms(q,1.4,-0.1,1.,0.2,0.8,0.3,0.4,3.,6.;np=64,nx=48)
        for z in (0.4+0.5im,3.0+0.8im,10.0+2im)
            v=testcausalgburegulatorchecks_RC.atom_value(g,z)
            @test v.contact_residual<1e-11
            @test abs(v.moment0)<1e-12
            @test v.value≈conj(testcausalgburegulatorchecks_RC.atom_value(rev,-conj(z)).value) atol=1e-12
        end
        if q==0
            original=testcausalgburegulatorchecks_RC.R.build_spectral_bubble(q,1.,0.2,1.4,-0.1,0.8;Phi=0.3,PhiBar=0.4,
                vacuum_cutoff_inv_fm=3.,thermal_cutoff_inv_fm=6.,momentum_nodes=64,angle_nodes=48)
            @test testcausalgburegulatorchecks_RC.atom_value(g,0.4+0.5im).value≈testcausalgburegulatorchecks_RC.R.spectral_bubble(original,0.1;eta_inv_fm=0.5).value atol=1e-12
            @test g.contact≈original.contact_inv_fm2 atol=1e-12
        end
        for l in (0.2,0.7,3.0,6.0)
            f(a,u,b,v,l,c)=testcausalgburegulatorchecks_RC.centered_cut(l,q,a,u,b,v,0.8,0.3,0.4,3.;component=c)
            @test f(1.,0.2,1.4,-0.1,l,:thermal)≈-f(1.4,-0.1,1.,0.2,-l,:thermal) atol=1e-12
            @test f(1.,0.,1.,0.,0.2,:vacuum)==0
        end
        @test abs(testcausalgburegulatorchecks_RC.centered_cut(0.3,q,1.,0.2,1.4,-0.1,0.8,0.3,0.4,6.;component=:thermal))<1e-12
    end
end
