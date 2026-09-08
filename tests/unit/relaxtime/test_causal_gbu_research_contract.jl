using Test
const _CGR_ROOT = normpath(joinpath(@__DIR__,"..","..",".."))
include(joinpath(_CGR_ROOT,"scripts","analysis","relaxtime","causal_gbu_research_utils.jl"))
const CGR = CausalGBUResearch

@testset "Frozen GBU research method" begin
    d = CGR.method_contract()
    @test !d["production_authorized"]
    @test d["observable"] == "fixed_quark_only_gbu_partial_yield"
    @test d["charge_to_baryon_ratio"] == 0.4
    @test d["strangeness_density_fm3"] == 0
    @test CGR.validate(CGR.Settings()) isa CGR.Settings
    @test_throws ArgumentError CGR.validate(CGR.Settings(lower=-0.1))
    @test_throws ArgumentError CGR.validate(CGR.Settings(nw=2))
    @test_throws ArgumentError CGR.validate(CGR.Settings(ne=2))
    for key in ("phase","weight","regulator","counting","background","bose_weight",
                "internal_frequency","q0_reference","observable","charge_to_baryon_ratio")
        changed=deepcopy(d)
        changed[key]="not_v1"
        @test_throws ErrorException CGR.validate_method(changed)
    end
    for key in ("phase_tolerance","root_residual_tolerance","complete_curve_requires_sparse_acceptance")
        changed=deepcopy(d)
        changed["gates"][key]=0
        @test_throws ErrorException CGR.validate_method(changed)
    end
    withenv("GBU_SPARSE_Q_NODES"=>"8","GBU_SPARSE_QMAX"=>"6.5","GBU_SPARSE_MESH"=>"64") do
        s=CGR.sparse_settings()
        @test s.nq==8 && s.qmax==6.5 && s.mesh==64
    end
    withenv("GBU_SPARSE_Q_NODES"=>"4") do
        @test_throws ArgumentError CGR.sparse_settings()
    end
    withenv("GBU_SPARSE_Q_NODES"=>"invalid") do
        @test_throws ArgumentError CGR.sparse_settings()
    end
    withenv("GBU_SPARSE_ENERGY_NODES"=>"128") do
        @test CGR.sparse_settings().ne==128
    end
    for shift in (-0.4,0.0,0.3), mass in (0.6,1.0,2.0), q in (0.0,0.5,1.0)
        w = hypot(mass,q)-shift
        @test CGR.q0_reference_coordinate(w,q,shift) ≈ mass-shift atol=1e-14
        # The Bose gap is k0, not k0-shift a second time.
        @test w ≈ hypot(CGR.q0_reference_coordinate(w,q,shift)+shift,q)-shift atol=1e-14
    end
    @test CGR.q0_reference_coordinate(0.2,1.0,0.3) === nothing
    @test_throws ArgumentError CGR.q0_reference_coordinate(1.0,-1.0,0.0)
end

include(joinpath(_CGR_ROOT,"scripts","analysis","relaxtime","audit_causal_gbu_routing.jl"))
@testset "Centered-regulator probe: q0 normalization and flavor reversal" begin
    P=CausalGBURouting
    for q in (0.0,0.5)
        g=P.centered_grid(q,1.0,0.2,1.4,-0.1,0.8,0.3,0.4,3.0,6.0;np=48,nx=32)
        reverse=P.centered_grid(q,1.4,-0.1,1.0,0.2,0.8,0.3,0.4,3.0,6.0;np=48,nx=32)
        @test P.value(g,0.6+0.4im) ≈ conj(P.value(reverse,-0.6+0.4im)) atol=1e-12
        if q==0
            old=CGR.build_spectral_bubble(q,1.0,0.2,1.4,-0.1,0.8;Phi=0.3,PhiBar=0.4,
                vacuum_cutoff_inv_fm=3.0,thermal_cutoff_inv_fm=6.0,momentum_nodes=48,angle_nodes=32)
            @test P.value(g,0.6+0.4im) ≈ CGR.spectral_bubble(old,0.3;eta_inv_fm=0.4).value atol=1e-12
            @test g.contact ≈ old.contact_inv_fm2 atol=1e-12
        end
    end
end

@testset "Independent contour count detects missed real roots" begin
    for roots in ([],[0.37],[0.37,0.37],[0.361,0.379],[-0.1,0.37,0.8,1.1])
        f(z) = prod((z-r for r in roots);init=1.0+0im)
        r = CGR.contour_count(f,0.0,1.0,-0.1,0.1)
        @test r.passed
        @test r.count == count(x->0<x<1,roots)
    end
    bad = CGR.contour_count(z->z,0.0,1.0,0.0,1.0)
    @test !bad.passed
    @test_throws ArgumentError CGR.contour_count(identity,1.0,0.0,-0.1,0.1)
    # Off-real-axis zeros are visible to the count but are not called bound states.
    unstable = CGR.contour_count(z->(z-0.5)^2+0.0025,0.0,1.0,-0.1,0.1)
    @test unstable.passed && unstable.count==2
    near_edge=CGR.contour_count(z->z-0.99999,0.0,1.0,-0.02,0.02)
    @test near_edge.passed && near_edge.count==1
    clustered=CGR.contour_count(z->z-(1-1e-8),0.0,1.0,-0.02,0.02)
    @test clustered.passed && clustered.count==1
    @test clustered.max_step<pi/2
    @test clustered.nodes_per_edge<=4096
end

@testset "Hard two-line domain has a linear small-q boundary term" begin
    m1,m2,L=1.0,1.4,3.0
    contact(q)=CGR.build_spectral_bubble(q,m1,0.0,m2,0.0,0.001;
        vacuum_cutoff_inv_fm=L,momentum_nodes=96,angle_nodes=48).contact_inv_fm2
    # The lost angular cap has integral int_{-1}^0 (-q*x) dx=q/2.
    expected=L^2*(1/hypot(L,m1)+1/hypot(L,m2))
    q=0.0001
    @test (contact(q)-contact(0.0))/q ≈ expected rtol=1e-4
    @test (contact(q/2)-contact(0.0))/(q/2) ≈ expected rtol=5e-5
    centered(q)=CausalGBURouting.centered_grid(q,m1,0.0,m2,0.0,0.001,0.0,0.0,L,L;
        np=96,nx=48).contact
    @test (centered(0.02)-centered(0.0))/(centered(0.01)-centered(0.0)) ≈ 4.0 rtol=1e-4
end

@testset "Threshold phase is a one-sided limit, not a fixed-offset count" begin
    threshold=2.0
    for sign in (-1,1)
        inverse(w)=sign-20im*sqrt(max(0.0,w-threshold))
        r=CGR.threshold_phase_limit(inverse,threshold)
        @test r.passed
        @test r.phase ≈ (sign<0 ? pi : 0.0) atol=3e-4
        @test r.offset==1e-10
        @test r.phase_change_over_pi<0.0005
    end
    @test !CGR.threshold_phase_limit(w->complex(w-threshold,0),threshold).passed
    @test !CGR.threshold_phase_limit(w->ComplexF64(NaN,0),threshold).passed
    @test !CGR.threshold_phase_limit(w->cis(0.5*sin(10log(max(w-threshold,1e-20)))),threshold).passed
    @test_throws ArgumentError CGR.threshold_phase_limit(identity,threshold;offsets=(1e-10,1e-7))
    @test_throws ArgumentError CGR.threshold_phase_limit(identity,1e10)
end
