using Test,LinearAlgebra
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","causal_gbu_source_continuation.jl"))
const SC=CausalGBUSourceContinuation

@testset "Fixed sea is not a Fermi-filled zero-temperature medium" begin
    @test SC.component_occupation(0.3,0.5,1,0.8,0.3,0.4,:vacuum)==0
    @test SC.S.occupation_component(-0.2,0.8,0.3,0.4,:vacuum)==1
    @test SC.reference_gap(0.,0.3,0.5,0.4,0.6)>0
    @test SC.reference_gap(0.,0.3,2.,0.4,0.)<0
    @test_throws ArgumentError SC.source_response(0.,0.3,2.,0.4,0.,0.8,2.,4.,[0.])
    @test_throws ArgumentError SC.source_response(0.,0.3,0.5,0.4,0.6,0.8,2.,4.,[0.];source_steps=[1.])
    for u in (0.1,0.6),ch in (:P,:S),component in (:vacuum,:thermal,:full)
        p=(0.,0.,0.1);r=(0.,0.,-0.1);a,b,v,T=0.3,0.4,u+0.04,0.8
        e1=hypot(0.1,a);e2=hypot(0.1,b)
        response=3sum(SC.S.spin_trace(p,r,a,b,s,t,ch)*SC.quotient(e1,u,s,e2,v,t,0.,T,0.3,0.4,component) for s in (-1,1),t in (-1,1))
        f(j)=SC.block_action(p,r,a,u,b,v,T,0.3,0.4,j,component;channel=ch)
        hs=(0.001,0.0005)
        ds=[(f(h)+f(-h)-2*f(0.))/h^2 for h in hs]
        @test abs(last(ds)+response)<1e-5
        @test abs(ds[2]-ds[1])<3e-5
        if u==0.1
            @test f(0.002)≈SC.S.block_action(p,r,a,u,b,v,T,0.3,0.4,0.002,component;channel=ch) atol=1e-13
        end
    end
end

@testset "Onset continuation preserves the existing polarization, not its regulator" begin
    for q in (0.,0.6),ch in (:P,:S)
        a,u,b,v,T,vc,tc=0.3,0.5,0.4,0.55,0.8,2.,4.
        zs=[0.,0.4+0.8im]
        s=SC.source_response(q,a,u,b,v,T,vc,tc,zs;Phi=0.3,PhiBar=0.4,nodes=48,channel=ch,source_steps=[0.002,0.001])
        g=SC.R.build_spectral_bubble(q,a,u,b,v,T;Phi=0.3,PhiBar=0.4,channel=ch,vacuum_cutoff_inv_fm=vc,thermal_cutoff_inv_fm=tc,momentum_nodes=96,angle_nodes=64)
        @test s.pi_inv_fm2[2]≈SC.R.spectral_bubble(g,0.4;eta_inv_fm=0.8).value atol=1e-7
        @test abs(last(s.source_curvatures_inv_fm2)+real(first(s.pi_inv_fm2)))<1e-6
        @test abs(diff(s.source_curvatures_inv_fm2)[1])<1e-6
        @test !s.below_individual_onset && !s.production_authorized
    end
end
