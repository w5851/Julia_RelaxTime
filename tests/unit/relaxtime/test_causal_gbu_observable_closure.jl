using Test,ForwardDiff
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","causal_gbu_observable_closure.jl"))
const OC=CausalGBUObservableClosure

@testset "GBU scalar selfenergy algebra permits signed ImPi" begin
    for K in (0.1,0.25),re in (-0.5,0.,0.6),impart in (-0.3,-1e-5,0.1),
        dp in (0.2+0.4im,-0.3+0.1im),branch in (-1,0,1)
        p=complex(re,impart)
        phase=-angle(1-4K*p)+2pi*branch
        r=OC.scalar_gbu_parts(p,dp,K;phase=phase)
        @test r.optical_error<1e-13
        @test r.weight_error<1e-13
        @test r.derivative_error<1e-13
        f(w)=OC.gbu_weight(-atan(imag(1-4K*(p+dp*w)),real(1-4K*(p+dp*w))))
        @test ForwardDiff.derivative(f,0.)≈r.derivative atol=1e-12
        @test !r.production_authorized
    end
    for d in (-1e-9,1e-8,1e-5)
        @test OC.gbu_weight(d)/(2d^3/3)≈1 rtol=1e-9
    end
    r=OC.scalar_gbu_parts(0.3+0.2im,0.,0.25)
    @test abs(r.weight-r.born_subtracted_phase)>1e-3
    @test abs(r.selfenergy_correction-r.full_product_imaginary)>1e-3
    @test_throws ArgumentError OC.scalar_gbu_parts(1.,0.,0.25)
    @test_throws ArgumentError OC.scalar_gbu_parts(0.2im,0.,0.25;phase=pi)
end

@testset "Nonlinear RPA must be applied after the same-source subtraction" begin
    r=OC.subtraction_logs(0.2+0.1im,0.3+0.2im,0.4+0.1im,0.25)
    @test abs(r.difference)>0.01
    @test !r.physical_equivalence_claimed
    # Vanishing medium correction is a special identity, not a general rule.
    r=OC.subtraction_logs(0.2+0.1im,0.4+0.1im,0.4+0.1im,0.25)
    @test abs(r.difference)<1e-14
end

@testset "Fixed-profile tagged count is not a total chemical derivative" begin
    edges=[0.03,0.3,1.,3.,8.]
    T=0.8
    phase(w,mu)=0.6w*exp(-w)*(1+0.7mu)
    delta(w)=phase(w,0.)
    terms=OC.partial_density_terms(delta,w->ForwardDiff.derivative(delta,w),
        w->0.7delta(w),edges,T)
    tag=ForwardDiff.derivative(nu->OC.tagged_pressure(delta,edges,T,nu),0.)
    total=ForwardDiff.derivative(mu->OC.tagged_pressure(w->phase(w,mu),edges,T,mu),0.)
    @test tag≈terms.tagged_density atol=1e-12
    @test total≈terms.total_chemical_derivative atol=1e-12
    # delta_mu>0 and sin(delta)^2*g>0 on this profile: the exact claim
    # is a positive omitted term, not an arbitrary minimum magnitude.
    @test terms.profile_response>0
    @test total-tag≈terms.profile_response atol=1e-12
    @test terms.integration_by_parts_error<1e-12
    @test abs(terms.endpoint_term)>1e-6
    @test !terms.fixed_profile_is_total_derivative
    # GBU does not turn a negative phase contribution into a positive count.
    neg=OC.partial_density_terms(w->-delta(w),w->-ForwardDiff.derivative(delta,w),
        w->0.,edges,T)
    @test neg.tagged_density<0
    # The simple-pole jump is preserved: W(n*pi+pi)-W(n*pi)=pi.
    for n in -2:2
        @test OC.gbu_weight((n+1)*pi)-OC.gbu_weight(n*pi)≈pi atol=2e-15
    end
end

@testset "GBU IR and small-phase tail estimates do not clip spectra" begin
    T,c=0.8,0.7
    for w in (1e-5,1e-6,1e-7)
        g=OC.bose(w,T,0.)
        kernel=g*(1+g)/T*OC.gbu_weight(c*w)
        @test kernel/w≈2T*c^3/3 rtol=1e-8
    end
    b=OC.tail_bound(0.8,0.03,5.,T,1.4)
    @test !b.rigorous_bound_certified
    @test !OC.tail_bound(0.8,0.03,5.,T,1.4;uniform_bounds_certified=true).rigorous_bound_certified
    @test OC.tail_bound(0.8,0.03,5.,T,1.4;uniform_bounds_certified=true,
        near_zero_branch_certified=true).rigorous_bound_certified
    for delta in range(-b.phase_bound,b.phase_bound;length=101)
        @test abs(OC.gbu_weight(delta))<=b.weight_bound
    end
    @test b.shell_bound_inv_fm2>0
    # The same inverse with a 2pi winding does NOT satisfy the small-tail bound.
    @test abs(OC.gbu_weight(2pi))>b.weight_bound
end

@testset "Signed Cauchy kernel has a conditional exterior zero-free bound" begin
    locations=[-2.,-0.4,1.7]; weights=[0.1,-0.3,0.2]
    S,M,K=maximum(abs,locations),sum(abs,weights),0.25
    radius=S+4K*M/pi+0.5
    b=OC.exterior_exclusion(S,M,K,radius;mass_bound_certified=true)
    @test b.exterior_zero_free_certified
    @test !OC.exterior_exclusion(S,M,K,radius).exterior_zero_free_certified
    @test !OC.exterior_exclusion(S,M,K,S+1e-3;mass_bound_certified=true).exterior_zero_free_certified
    for theta in range(0,pi;length=65)
        z=radius*cis(theta)
        p=sum(w/(x-z) for (x,w) in zip(locations,weights))/pi
        @test abs(p)<=b.bubble_bound_inv_fm2
        @test abs(1-4K*p)>=b.inverse_lower_bound
    end
    @test_throws ArgumentError OC.exterior_exclusion(S,M,K,S)
end

@testset "Signed auxiliary Gaussian determinant: contour is conditional" begin
    for E in (0.8,2.1),P in (1.2,1.8),T in (0.6,0.9)
        r=OC.toy_matsubara(E,P,T)
        @test abs(r.matsubara-r.exact)<r.remainder_bound+1e-13
        @test !r.is_project_result
        @test sign(r.signed_excess_count)==sign(P-E)
        f(z)=OC.toy_inverse(z,E^2,P^2)
        count=OC.R.contour_count(f,-3.,3.,0.01,3.;nodes=64)
        @test count.passed && count.count==0
    end
    unstable(z)=OC.toy_inverse(z,-0.36,1.44)
    count=OC.R.contour_count(unstable,-3.,3.,0.01,3.;nodes=64)
    @test count.passed && count.count==1
    @test unstable(0.)<0
    @test_throws ArgumentError OC.toy_matsubara(-0.6,1.2,0.8)
end
