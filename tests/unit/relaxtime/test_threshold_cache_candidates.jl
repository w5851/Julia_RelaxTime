using Test

if !isdefined(Main,:ThresholdCacheCandidates)
    include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","ThresholdCacheCandidates.jl"))
end
using .ThresholdCacheCandidates: build_adaptive_sigma, composite_sigma, PairLoopPlan,
    pair_b0_pv, functional_rate, paired_cross_section
using .ThresholdUpgradeCommon

@testset "Weighted interpolation resolves an integrable threshold" begin
    exact=build_adaptive_sigma(x->1/x,[1e-8,0.6];representation=:x_sigma,rtol=1e-8)
    @test exact.converged
    @test length(exact.nodes)==2
    @test exact.estimated_error<1e-14
    @test isapprox(exact(0.02),50.0;rtol=1e-14)
    @test_throws DomainError exact(0.0)
    @test_throws DomainError exact(0.7)
    @test_throws ArgumentError build_adaptive_sigma(identity,[0.0,1.0])
    @test_throws ArgumentError build_adaptive_sigma(identity,[1.0,2.0];rtol=NaN)
    @test_throws DomainError build_adaptive_sigma(x->-1.0,[0.1,1.0])

    a=1e-4; lo=1e-9; hi=0.6
    sigma(x)=1/(a+x)+0.02
    seeds=[lo,a/4,a,4a,0.01,hi]
    truth=hi-lo-a*log((a+hi)/(a+lo))+0.01*(hi^2-lo^2)
    weighted=build_adaptive_sigma(sigma,seeds;rtol=1e-4)
    raw=build_adaptive_sigma(sigma,seeds;representation=:raw,rtol=1e-4)
    @test weighted.converged && raw.converged
    integral=sum((weighted.nodes[i+1]-weighted.nodes[i])*
        (weighted.values[i+1]+weighted.values[i])/2 for i in 1:length(weighted.nodes)-1)
    @test abs(integral/truth-1)<1.1e-4
    @test length(weighted.nodes)<length(raw.nodes)
    @test all(isapprox(weighted(x),sigma(x);rtol=1e-14) for x in weighted.nodes)

    exhausted=build_adaptive_sigma(sigma,[lo,hi];representation=:raw,rtol=1e-12,max_nodes=2)
    @test !exhausted.converged
    @test_throws ArgumentError composite_sigma(exhausted,1.0,identity,identity)
end

@testset "Paired PV preserves finite threshold and square-root cut" begin
    vacuum=(T=0.001,Φ=0.0,Φbar=0.0,ξ=0.0)
    plan=PairLoopPlan(1.0,0.0,vacuum;n_nodes=128)
    @test isapprox(real(pair_b0_pv(plan,4.0)),2asinh(plan.cutoff);atol=2e-13)
    @test imag(pair_b0_pv(plan,4.0))==0
    @test_throws DomainError pair_b0_pv(plan,4.0+4plan.cutoff^2)
    @test_throws ArgumentError pair_b0_pv(plan,-1.0)
    @test_throws ArgumentError PairLoopPlan(0.0,0.0,vacuum)

    tp=(T=0.8,Φ=0.5,Φbar=0.5,ξ=0.3)
    p64=PairLoopPlan(1.0,0.15,tp;n_nodes=64)
    p128=PairLoopPlan(1.0,0.15,tp;n_nodes=128)
    for x in (-0.4,-1e-4,-1e-8,0.0,1e-8,1e-4,0.4)
        @test isapprox(pair_b0_pv(p64,4+x),pair_b0_pv(p128,4+x);rtol=2e-10,atol=1e-12)
    end
    x=1e-8
    slope_below=(real(pair_b0_pv(p128,4-x))-real(pair_b0_pv(p128,4.0)))/sqrt(x)
    slope_above=imag(pair_b0_pv(p128,4+x))/sqrt(x)
    @test isapprox(slope_below,-π*p128.F0/2;rtol=1e-3)
    @test isapprox(slope_above,π*p128.F0/2;rtol=1e-6)
    # Force a quadrature-node/pole coincidence to exercise the removable limit.
    @test isfinite(pair_b0_pv(p64,4+4p64.p2[20]))
end

@testset "Callable thermal consumer and temporary loop injection" begin
    row=(m_u=1.4,m_d=1.4,m_s=2.2,muq_fm=0.0,T_fm=0.8,Phi=0.5,Phibar=0.5,xi=0.2)
    state=retained_state(row)
    qp,tp,kc=state.quark_params,state.thermo_params,state.K_coeffs
    problem=AS._prepare_thermal_cache_problem(:uubar_to_uubar,qp,tp,kc;
        p_nodes=8,angle_nodes=4,phi_nodes=6,sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm)
    measure=AS.ThermalWeightedSampling.thermal_measure(problem.plan,problem.context.bounds)
    reference=0.75*sum(measure.weights)
    candidate=functional_rate(s->0.75,:uubar_to_uubar,state;p_nodes=8,angle_nodes=4,phi_nodes=6)
    @test isapprox(candidate,reference;rtol=1e-13,atol=0)

    s=4qp.m.u^2+0.1
    baseline=CS.total_cross_section(:uubar_to_uubar,s,qp,tp,kc;n_points=4)
    plans=(u=PairLoopPlan(qp.m.u,qp.μ.u,tp),s=PairLoopPlan(qp.m.s,qp.μ.s,tp))
    changed=paired_cross_section(:uubar_to_uubar,s,state,plans;n_points=4)
    @test isfinite(changed) && changed>=0
    @test CS.total_cross_section(:uubar_to_uubar,s,qp,tp,kc;n_points=4)==baseline
    other=merge(state,(thermo_params=merge(tp,(ξ=0.3,)),))
    @test_throws ArgumentError paired_cross_section(:uubar_to_uubar,s,other,plans;n_points=4)
end
