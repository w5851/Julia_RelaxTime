using Test
const _HEALTH_ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
include(joinpath(_HEALTH_ROOT,"scripts","analysis","relaxtime","causal_gbu_pnjl_stability.jl"))
const H=CausalGBUPNJLStability
include(joinpath(_HEALTH_ROOT,"scripts","analysis","relaxtime","summarize_causal_gbu_pnjl_stability.jl"))
const HS=CausalGBUPNJLStabilitySummary

@testset "Health evidence rejects contour promotion and false passes" begin
    row=(coarse_winding=0.,fine_winding=0.,max_step=0.2,minimum_inverse_abs=0.1,
        count=0,resolved=true,full_UHP_certified=false,unresolved_near_axis_strip=true)
    @test HS.check_contour(row)
    @test_throws ErrorException HS.check_contour(merge(row,(max_step=1.,)))
    @test_throws ErrorException HS.check_contour(merge(row,(full_UHP_certified=true,)))
end

@testset "PNJL occupation consistency and gap algebra" begin
    state=H.tadpole_gap_state((1.,2.,3.),ntuple(_->-4pi^2/3,3),(.1,.2,.3),.2,.1)
    @test collect(state.condensates_inv_fm3)≈[-1.,-2.,-3.]
    @test collect(state.masses_inv_fm)≈[2.1,2.4,3.1]
    @test state.K12_fm2≈0.35 && state.K45_fm2≈0.3
    for phi in (0.,0.3,1.),bar in (0.,0.4,1.),x in (-100.,-2.,0.,2.,100.)
        @test H.P.occupation(x,phi,bar)+H.P.occupation(-x,bar,phi)≈1 atol=3e-16
        @test H.P.occupation(x+0.01,phi,bar)<=H.P.occupation(x,phi,bar)
    end
    for phi in (0.,0.3,1.),mu in (-0.2,0.2)
        m,T,L=1.3,0.8,3.
        ps,ws=H.R.gauleg(0.,16.,160)
        A=Main.RelaxTime.OneLoopIntegrals.A(m,mu,T,phi,0.4,ps,ws)
        g=H.R.build_spectral_bubble(0.,m,mu,m,mu,T;Phi=phi,PhiBar=0.4,
            vacuum_cutoff_inv_fm=Main.Constants_PNJL.Λ_inv_fm,thermal_cutoff_inv_fm=16.,
            momentum_nodes=160,angle_nodes=16)
        pi0=sum(g.pi_residues_inv_fm3 ./ (-g.poles_inv_fm))
        @test g.contact_inv_fm2≈2A atol=2e-11
        @test pi0≈-3A/(4pi^2) atol=2e-11
        K=-pi^2/(3A)
        @test abs(1-4K*pi0)<1e-11
    end
end

# Loading the runner checks executable syntax without launching diagnostics.
include(joinpath(_HEALTH_ROOT,"scripts","analysis","relaxtime","audit_causal_gbu_pnjl_stability.jl"))
@testset "PNJL health runner settings" begin
    withenv("GBU_HEALTH_Q"=>"0,1.4","GBU_HEALTH_MESH"=>"128","GBU_HEALTH_NODES"=>"256") do
        s=CausalGBUPNJLStabilityAudit.settings()
        @test s.qs==[0.,1.4] && s.mesh==128 && s.np==256
    end
    withenv("GBU_HEALTH_Q"=>"1.4") do
        @test_throws ArgumentError CausalGBUPNJLStabilityAudit.settings()
    end
    withenv("GBU_HEALTH_CHANNELS"=>"unknown") do
        @test_throws ArgumentError CausalGBUPNJLStabilityAudit.settings()
    end
    withenv("GBU_HEALTH_VARIANTS"=>"thermal24","GBU_HEALTH_CHANNELS"=>"pi_plus") do
        s=CausalGBUPNJLStabilityAudit.settings()
        @test s.channels==[:pi_plus] && s.variants==["thermal24"]
    end
end

@testset "Spectral derivative, signed residues and passivity" begin
    p=H.R.PiecewiseSpectralFunction([-3.,-2.,-1.,0.,1.,2.,3.],[0.,-1.,0.,0.,0.,1.,0.])
    for z in (0.,0.3,4.,0.4+0.6im,1e5+1im)
        h=z==1e5+1im ? 0.1 : 1e-5
        fd=(H.R.cauchy_transform(p,z+h)-H.R.cauchy_transform(p,z-h))/(2h)
        @test H.spectral_derivative(p,z)≈fd rtol=2e-8 atol=2e-10
    end
    @test_throws ArgumentError H.spectral_derivative(p,2.)
    @test H.passivity_audit(p,0.01,0.).sign_sufficient_condition_float64
    bad=H.R.PiecewiseSpectralFunction([-2.,-1.,0.,1.,2.],[0.,-1.,0.,-1.,0.])
    @test !H.passivity_audit(bad,0.01,0.).numerical_passivity
    # The exact zero endpoints must stay zero during the sign check.
    tail=H.R.PiecewiseSpectralFunction([0.,3.12345,3.23456,4.98765,6.],
        [0.,0.,1.34567,0.,0.])
    @test H.passivity_audit(tail,0.001,0.).minimum_k0_rho==0
    # An interior negative product cannot be hidden by positive endpoints.
    crossing=H.R.PiecewiseSpectralFunction([-2.,-1.,1.,2.],[0.,-2.,1.,0.])
    @test H.passivity_audit(crossing,0.001,0.).minimum_k0_rho<0
    @test H.pole_weight(-2.,0.2,1.).sign_passed
    @test H.pole_weight(2.,0.2,-1.).sign_passed
    @test !H.pole_weight(2.,0.2,1.).sign_passed
    @test !H.pole_weight(0.,0.2,1.).simple
    near=z->log(1-z)
    adaptive=H.gap_derivative_check(near,0.9999,-1.,1.,-1/(1-0.9999))
    @test adaptive.passed && adaptive.error<1e-5 && adaptive.halvings>1
    @test adaptive.initial_error>1e-5 && adaptive.step_change<1e-5
    limited=H.gap_derivative_check(near,0.9999,-1.,1.,-1/(1-0.9999);max_halvings=1)
    @test !limited.passed
    @test_throws ArgumentError H.gap_derivative_check(near,1.,-1.,1.,0.)
end

@testset "UHP contour detects instabilities without counting real poles" begin
    for alpha in (0.2,0.7)
        f(z)=(z^2+alpha^2)/(z^2-1)
        r=H.upper_rectangle_count(f,3.,0.01)
        @test r.resolved && r.count==1
        @test r.unresolved_near_axis_strip && !r.full_UHP_certified
    end
    stable=z->(z^2-0.25)/(z^2-1)
    @test H.upper_rectangle_count(stable,3.,0.01).count==0
    realroot=0.21713
    near=z->(z-realroot)*(1+0.2z)
    anchored=H.upper_rectangle_count(near,3.,0.003;
        anchors=[realroot-0.003,realroot,realroot+0.003],nodes=48,max_nodes=768)
    @test anchored.resolved && anchored.count==0 && anchored.max_step<pi/4
    @test anchored.contour_nodes>2*(4*48+3)
    @test_throws ArgumentError H.upper_rectangle_count(stable,3.,0.)
end

@testset "Single-sphere exact domain, routing and q0 parity" begin
    zs=[0.4+0.6im,3.2+0.8im]
    args=(1.,0.2,1.4,-0.1,0.8,0.3,0.4,3.,8.,zs)
    for q in (0.,0.8,3.2)
        a=H.single_sphere_loop(q,args...;nodes=80,coordinates=:angular)
        b=H.single_sphere_loop(q,args...;nodes=80,coordinates=:radial)
        @test a.values≈b.values atol=2e-9
        @test a.contact≈b.contact atol=2e-10
        @test !a.production_authorized
    end
    zero=H.single_sphere_loop(0.,args...;nodes=96)
    one=H.single_sphere_loop(0.8,args...;nodes=96)
    @test one.A1≈zero.A1 atol=2e-10
    @test abs(one.A2-zero.A2)>0.01 # even without the third condition A2 is q-dependent
    reverse=H.single_sphere_loop(0.8,1.4,-0.1,1.,0.2,0.8,0.3,0.4,3.,8.,-conj.(zs);nodes=96)
    @test maximum(abs.(one.values.-conj.(reverse.values)))>1e-4
    @test_throws ArgumentError H.single_sphere_loop(0.,args...;coordinates=:unknown)
end
