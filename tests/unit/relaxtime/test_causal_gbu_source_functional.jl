using Test,LinearAlgebra,ForwardDiff
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","causal_gbu_source_functional.jl"))
const testcausalgbusourcefunctional_SF=CausalGBUSourceFunctional

@testset "Fixed projector determinant differentiates to two projected lines" begin
    d=[1.,2.,3.,4.]
    V=[0. 1. 2. 0.;1. 0. 0. 3.;2. 0. 0. 4.;0. 3. 4. 0.]
    for keep in (falses(4),trues(4),Bool[1,1,0,0],Bool[1,0,1,0])
        f=j->testcausalgbusourcefunctional_SF.projected_action(d,V,keep,j)
        exact=testcausalgbusourcefunctional_SF.projected_hessian(d,V,keep)
        @test ForwardDiff.derivative(j->ForwardDiff.derivative(f,j),0.)≈exact atol=1e-13
        @test exact>=0
    end
    keep=Bool[1,1,0,0]
    two=testcausalgbusourcefunctional_SF.projected_hessian(d,V,keep)
    naive_one=sum(abs2(V[i,j])/(d[i]*d[j]) for i in 1:4,j in 1:4 if keep[i])
    @test naive_one>two
    # Same-field subtraction: the large vacuum determinant cancels exactly.
    f=j->testcausalgbusourcefunctional_SF.projected_action(d,V,keep,j)+testcausalgbusourcefunctional_SF.projected_action(d.+0.7,V,trues(4),j)-
        testcausalgbusourcefunctional_SF.projected_action(d,V,trues(4),j)
    @test ForwardDiff.derivative(j->ForwardDiff.derivative(f,j),0.)≈
        two+testcausalgbusourcefunctional_SF.projected_hessian(d.+0.7,V,trues(4))-testcausalgbusourcefunctional_SF.projected_hessian(d,V,trues(4)) atol=1e-13
    @test_throws ArgumentError testcausalgbusourcefunctional_SF.projected_action([-1.],ones(1,1),[true],0.)
end

@testset "Complex color determinant reproduces PNJL algebra, not positivity" begin
    for (phi,bar) in ((0.,0.),(0.3,0.4),(1.,1.)),x in (-2.,0.,3.)
        ell=testcausalgbusourcefunctional_SF.color_eigenvalues(phi,bar)
        y=exp(-x)
        @test prod(1 .+ell*y)≈1+3phi*y+3bar*y^2+y^3 rtol=2e-13
        @test sum(ell*y./(1 .+ell*y))/3≈testcausalgbusourcefunctional_SF.P.occupation(x,phi,bar) atol=1e-13
        @test prod(1 .+y./ell)≈1+3bar*y+3phi*y^2+y^3 rtol=2e-13
    end
    @test any(abs.(abs.(testcausalgbusourcefunctional_SF.color_eigenvalues(0.3,0.4)).-1).>0.01)
    @test_throws ArgumentError testcausalgbusourcefunctional_SF.color_eigenvalues(-0.1,0.3)
end

@testset "Dirac source trace equals the on-shell A/B numerator" begin
    @test testcausalgbusourcefunctional_SF.VP≈testcausalgbusourcefunctional_SF.VP'
    @test testcausalgbusourcefunctional_SF.VP^2≈testcausalgbusourcefunctional_SF.ID4
    for p in ((0.,0.,0.),(0.3,0.2,0.7)),r in ((0.2,0.,-0.1),(0.,0.,0.)),
        (a,b) in ((1.2,1.2),(1.1,1.9)),s in (-1,1),t in (-1,1),ch in (:P,:S)
        e1=sqrt(sum(abs2,p)+a^2); e2=sqrt(sum(abs2,r)+b^2)
        q2=sum((p[i]-r[i])^2 for i in 1:3)
        m2=ch===:P ? (a-b)^2 : (a+b)^2
        formula=-s*t*((s*e1-t*e2)^2-q2-m2)/(2*e1*e2)
        @test testcausalgbusourcefunctional_SF.spin_trace(p,r,a,b,s,t,ch)≈formula atol=2e-15
        @test testcausalgbusourcefunctional_SF.spin_trace(p,r,a,b,s,t,ch)>=-1e-15
    end
end

@testset "Static source trace-log and ordered dynamic Hessian" begin
    for ch in (:P,:S),component in (:vacuum,:thermal),phi in (0.3,1.)
        p=(0.2,0.,0.5);r=(0.2,0.,-0.4)
        a,u,b,v,T=1.2,0.1,1.7,0.25,0.8
        e1=sqrt(sum(abs2,p)+a^2);e2=sqrt(sum(abs2,r)+b^2)
        direct=3sum(testcausalgbusourcefunctional_SF.spin_trace(p,r,a,b,s,t,ch)*
            testcausalgbusourcefunctional_SF.response_quotient(s*e1-u,t*e2-v,0.,T,phi,0.4,component)
            for s in (-1,1),t in (-1,1))
        f(j)=testcausalgbusourcefunctional_SF.block_action(p,r,a,u,b,v,T,phi,0.4,j,component;channel=ch)
        h=0.002
        d1=(f(h)+f(-h)-2f(0.))/h^2
        d2=(f(h/2)+f(-h/2)-2f(0.))/(h/2)^2
        @test abs(d2+direct)<2e-6
        @test abs(d2-d1)<2e-6
    end
    for phi in (0.3,1.),tc in (2.5,6.)
        a,u,T=1.2,0.1,0.8
        response=testcausalgbusourcefunctional_SF.source_response(0.,a,u,a,u,T,2.5,tc,[0.];Phi=phi,PhiBar=0.4,nodes=32)
        state=testcausalgbusourcefunctional_SF.B.static_flavor(a,u,T,phi,0.4,2.5,tc;nodes=96)
        @test response.pi_inv_fm2[1]≈-state.condensate_inv_fm3/a atol=1e-9
        @test !response.positive_physical_trace_certified && !response.production_authorized
    end
    @test_throws ArgumentError testcausalgbusourcefunctional_SF.source_response(0.,1.,1.1,1.,0.,0.8,2.,4.,[0.])
    @test_throws ArgumentError testcausalgbusourcefunctional_SF.source_response(0.,1.,0.,1.,0.,0.8,2.,4.,[1.])
    # Exercise the integrated finite-difference assembly, not only a local block.
    for q in (0.,0.6),ch in (:P,:S)
        rs=testcausalgbusourcefunctional_SF.source_response(q,1.2,0.1,1.7,0.25,0.8,2.5,4.,[0.,1+0.8im];
            Phi=0.3,PhiBar=0.4,nodes=24,channel=ch,source_steps=[0.004,0.002,0.001])
        @test abs(rs.source_curvatures_inv_fm2[end]+real(rs.pi_inv_fm2[1]))<1e-6
        @test abs(rs.source_curvatures_inv_fm2[end]-rs.source_curvatures_inv_fm2[end-1])<1e-6
        g=testcausalgbusourcefunctional_SF.R.build_spectral_bubble(q,1.2,0.1,1.7,0.25,0.8;Phi=0.3,PhiBar=0.4,
            channel=ch,vacuum_cutoff_inv_fm=2.5,thermal_cutoff_inv_fm=4.,
            momentum_nodes=48,angle_nodes=48)
        @test testcausalgbusourcefunctional_SF.R.spectral_bubble(g,1.;eta_inv_fm=0.8).value≈rs.pi_inv_fm2[2] atol=1e-8
    end
end
