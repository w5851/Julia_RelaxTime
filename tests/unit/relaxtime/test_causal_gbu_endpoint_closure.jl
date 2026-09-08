using Test
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","causal_gbu_endpoint_closure.jl"))
const testcausalgbuendpointclosure_EC=CausalGBUEndpointClosure

@testset "Projected exterior bound independent of spectral interpolation" begin
    @test testcausalgbuendpointclosure_EC.lens_volume(0.,3.)≈4pi*3^3/3
    @test testcausalgbuendpointclosure_EC.lens_volume(6.,3.)==0
    @test testcausalgbuendpointclosure_EC.lens_volume(7.,3.)==0
    for q in (0.,0.5,3.,7.,21.)
        r=testcausalgbuendpointclosure_EC.projected_exterior(q,1.,2.,3.,10.,0.25)
        @test r.inverse_deviation_bound<=0.5+1e-14
        @test r.analytic_exterior_exclusion
        @test !r.directed_rounding_certified
    end
end

@testset "Signed endpoint toy: separate root, continuum and weak eta limit" begin
    a,S,c,K,T=1.,2.,1.,0.7,0.8
    root=testcausalgbuendpointclosure_EC.toy_exterior_root(a,S,c,K).root
    F(z)=testcausalgbuendpointclosure_EC.toy_inverse(z,a,S,c,K)
    Fp(z)=4K*c/pi*(-1/(S-z)+1/(S+z)+1/(a-z)-1/(a+z))
    phase(w)=-angle(F(complex(w)))
    W(w)=testcausalgbuendpointclosure_EC.O.gbu_weight(phase(w))
    dW(w)=2sin(phase(w))^2*(-imag(Fp(complex(w))/F(complex(w))))
    g(w)=testcausalgbuendpointclosure_EC.O.bose(w,T,0.)
    beta(w)=g(w)*(1+g(w))/T
    integrate(f,es)=testcausalgbuendpointclosure_EC.A.integrate_panels(f,es;nodes=16,atol=1e-9,max_depth=24)
    tag=integrate(w->beta(w)*W(w)/pi,[a,S,root])
    # Subtract the line joining g at both endpoints, retaining its exact terms.
    line(w)=g(a)+(g(S)-g(a))*(w-a)/(S-a)
    slope=(g(S)-g(a))/(S-a)
    residual=integrate(w->((g(w)-line(w))*dW(w)-slope*W(w))/pi,[a,S])
    continuum=residual.value-g(S)
    discrete=g(root)
    @test tag.converged && residual.converged
    @test tag.value≈continuum+discrete atol=1e-8
    @test tag.value<0 && continuum<0 && discrete>0
    # An auxiliary root cannot be omitted even though its residue is negative.
    @test abs(tag.value-continuum)>0.01
    etas=(0.02,0.005,0.00125)
    differences=Float64[]
    for eta in etas
        edges=sort!(unique!([0.2,a-eta,a,a+eta,S-eta,S,S+eta,root-eta,root,root+eta,4.]))
        smooth=integrate(w->beta(w)*testcausalgbuendpointclosure_EC.O.gbu_weight(-angle(F(w+eta*im)))/pi,edges)
        @test smooth.converged
        push!(differences,abs(smooth.value-tag.value))
    end
    @test differences[3]<differences[2]<differences[1]
    @test differences[3]<0.001
    # This test concerns weighted integrals, not equality at poles/cut endpoints.
end

@testset "Finite thermal edge forces exterior zeros, not UHP instability" begin
    for phi in (0.,0.3,1.),bar in (0.,0.4,1.),T in (0.6,0.9),L in (10.,20.,24.)
        ep=testcausalgbuendpointclosure_EC.thermal_endpoint(1.6,0.3,2.6,0.01,T,phi,bar,3.,L)
        @test ep.rho_inside_inv_fm2<0
        r=testcausalgbuendpointclosure_EC.endpoint_obstruction(ep.rho_inside_inv_fm2,ep.lambda_endpoint_inv_fm,1.,-0.1,0.24,T,0.29)
        @test r.at_least_one_positive_exterior_zero
        @test r.sub_float64_resolution
        @test !r.UHP_instability_implied
        @test !r.bose_bound_is_total_density_bound
        @test !r.uniqueness_proved
    end
    @test_throws ArgumentError testcausalgbuendpointclosure_EC.thermal_endpoint(1.,0.,2.,0.,0.8,1.,1.,3.,3.)
    @test_throws ArgumentError testcausalgbuendpointclosure_EC.endpoint_obstruction(0.1,5.,1.,0.,0.25,0.8,0.)
end

@testset "Exact signed-band counterexample has a real negative-residue root" begin
    for c in (0.3,0.6,1.),K in (0.4,0.7)
        r=testcausalgbuendpointclosure_EC.toy_exterior_root(1.,2.,c,K)
        @test r.root>2
        @test abs(testcausalgbuendpointclosure_EC.toy_inverse(r.root,1.,2.,c,K))<1e-10
        @test r.residue<0
        @test !r.is_positive_physical_bound_state
        @test r.phase_jump≈pi
        @test real(testcausalgbuendpointclosure_EC.toy_inverse(1im,1.,2.,c,K))>1
        step=min(1e-4,(r.root-2)/10)
        before=-angle(testcausalgbuendpointclosure_EC.toy_inverse(r.root-step+step*1e-6im,1.,2.,c,K))
        after=-angle(testcausalgbuendpointclosure_EC.toy_inverse(r.root+step+step*1e-6im,1.,2.,c,K))
        @test after-before≈pi atol=1e-5
        # The continuum ends at -pi; the exterior zero supplies +pi.
        # Counting only positive-residue physical states would miss this unit.
        continuum_units=(testcausalgbuendpointclosure_EC.O.gbu_weight(-pi)-testcausalgbuendpointclosure_EC.O.gbu_weight(0.))/pi
        @test continuum_units+1≈0 atol=1e-14
        @test continuum_units+0≈-1 atol=1e-14
    end
    s=testcausalgbuendpointclosure_EC.stage_status(analytic_structure=false,physical_counting=true,integral_convergence=true)
    @test s.step2=="blocked_by_step1" && s.step3=="blocked_by_counting"
    @test !s.research_production_ready && !s.production_authorized
end
