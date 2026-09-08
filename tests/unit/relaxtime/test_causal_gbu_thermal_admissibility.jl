using Test
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","causal_gbu_thermal_admissibility.jl"))
const A=CausalGBUThermalAdmissibility

@testset "PNJL positive polynomial and thermal-only pair signs" begin
    for phi in (0.,0.3,1.),bar in (0.,0.4,1.),x in (-20.,-2.,0.,2.,20.)
        s=A.occupation_statistics(x,phi,bar)
        @test s.occupation≈A.P.occupation(x,phi,bar) atol=3e-16
        @test s.derivative<=0
        h=1e-5
        @test s.derivative≈(A.P.occupation(x+h,phi,bar)-A.P.occupation(x-h,phi,bar))/(2h) atol=1e-9
        for y in (-3.,1.,4.)
            pair=1-A.P.occupation(x,phi,bar)-A.P.occupation(y,bar,phi)
            @test (x+y)*pair>=-1e-14
        end
    end
    for phi in (0.,0.3,1.),p in (3.01,3.5,5.9)
        r=A.q0_pair(p,1.,0.2,1.4,-0.1,0.8,phi,0.4,3.,6.)
        @test r.thermal_only_interval && r.k0_inv_fm>0 && r.imaginary_inv_fm2<0
        hard=A.q0_pair(p,1.,0.2,1.4,-0.1,0.8,phi,0.4,3.,3.)
        @test hard.imaginary_inv_fm2==0
    end
    @test_throws ArgumentError A.q0_pair(4.,1.,0.2,1.4,-0.1,0.8,0.3,0.4,3.,2.)
end

@testset "Partition hints and runner scope are explicit" begin
    for q in (0.,1.2,6.2)
        g=A.R.build_spectral_bubble(q,1.,0.2,1.4,-0.1,0.8;Phi=0.3,PhiBar=0.4,
            vacuum_cutoff_inv_fm=3.,thermal_cutoff_inv_fm=6.,momentum_nodes=8,angle_nodes=8)
        edges=A.cut_panels(g)
        @test all(diff(edges).>0)
        @test first(edges)==-last(edges)==-(hypot(1.,6.)+hypot(1.4,6.))
        @test 0. in edges && 0.2-(-0.1) in edges
    end
    include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime",
        "audit_causal_gbu_thermal_admissibility.jl"))
    s=CausalGBUThermalAdmissibilityAudit.settings()
    @test s.channels==A.R.CHANNELS
    @test s.occupations==("fermi","pnjl")
    @test 0. in s.witness_q && 0. in s.etas
    @test minimum(filter(>(0),s.etas))==0.0003
    @test s.direct_node_tolerance_inv_fm2==1e-6
    include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime",
        "summarize_causal_gbu_thermal_admissibility.jl"))
    w=(k0_inv_fm=2.,vacuum_pair_inv_fm2=0.,landau_inv_fm2=0.,thermal_pair_inv_fm2=-0.2,
        original_cut_inv_fm2=-0.2,hard_cut_inv_fm2=0.,cut_error_inv_fm2=1e-12,
        node_change_inv_fm2=1e-12,q_inv_fm=0.,analytic_q0_error_inv_fm2=1e-12)
    v=CausalGBUThermalAdmissibilitySummary.witness_verdict
    @test v(w,1e-9)==(negative=true,numeric=true)
    @test !v(merge(w,(cut_error_inv_fm2=NaN,)),1e-9).numeric
    @test !v(merge(w,(original_cut_inv_fm2=0.2,)),1e-9).negative
end

@testset "Continuous Cauchy oracle resolves PV and near-axis kernels" begin
    # rho(v)=v on [-1,1], with supplied discontinuities at the outer endpoints.
    rho=x->x
    for z in (0.,0.2,0.2+1e-2im,0.2+1e-5im,2.0+0.1im)
        logpart=imag(z)==0 ? complex(log(abs((1-z)/(-1-z))),pi) :
            log(1-z)-log(-1-z)
        exact=(2+z*logpart)/pi
        r=A.continuous_cauchy(rho,[-1.,0.,1.],z;atol=1e-9)
        @test r.value≈exact atol=3e-9
        @test r.converged && !r.rigorous_error_bound
    end
    q=A.integrate_panels(x->sqrt(x),[0.,1.];atol=1e-8)
    @test q.value≈2/3 atol=1e-8
    failed=A.integrate_panels(x->sqrt(x),[0.,1.];nodes=4,atol=1e-15,max_depth=0)
    @test !failed.converged && !failed.rigorous_error_bound
    # A narrow panel must not be forced to a vanishing fraction of the GLOBAL
    # error budget. Its observed error is still counted without clipping.
    localized=x->x<1e-10 ? sqrt(x/1e-10) : 0.
    narrow=A.integrate_panels(localized,[0.,1e-10,1.];nodes=4,atol=1e-9,max_depth=0)
    @test narrow.converged && narrow.error_estimate>0
    @test narrow.value≈2e-10/3 atol=1e-13
    @test_throws ArgumentError A.continuous_cauchy(rho,[-1.,1.],0.2-1im)
    @test_throws ArgumentError A.continuous_cauchy(rho,[1.,-1.],0.2)
    @test_throws ArgumentError A.continuous_cauchy(rho,[-1.,1.],1.)
    # An INTERNAL hard-cutoff jump is integrated with one-sided panels.
    jump=x->x<0 ? 1. : -0.2
    for z in (0.4,0.4+1e-4im)
        lp(a,b)=imag(z)==0 ? complex(log(abs((b-z)/(a-z))),a<real(z)<b ? pi : 0.) :
            log(b-z)-log(a-z)
        expected=(lp(-1.,0.)-0.2lp(0.,1.))/pi
        actual=A.continuous_cauchy(jump,[-1.,0.,1.],z;atol=1e-9)
        @test actual.value≈expected atol=3e-9
        @test actual.converged
    end
end

@testset "Signed RPA spectrum is not removed by a real contact shift" begin
    for re in (-2.,0.,2.),impart in (-0.1,0.1)
        r=A.rpa_spectral(complex(re,impart),0.2)
        @test r.imaginary≈r.identity_value atol=1e-15
        @test sign(r.imaginary)==sign(impart)
    end
    @test_throws ArgumentError A.rpa_spectral(1.25,0.2)
    @test !A.count_transfer_margin(1.,1e-5).count_transfer_certified
    @test A.count_transfer_margin(1.,1e-5;uniform_bound_certified=true).count_transfer_certified
    @test !A.count_transfer_margin(1.,2.;uniform_bound_certified=true).count_transfer_certified
end
