using Test
if !isdefined(Main,:ThermalAverageCandidates)
    include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","ThermalAverageCandidates.jl"))
end
using .ThermalAverageCandidates
using .ThresholdUpgradeCommon

@testset "Exact thermal measure across all channels and asymmetric grids" begin
    ts=AS.ThermalWeightedSampling
    for xi in (-0.48,0.0,0.49), custom in (false,true)
        qp=(m=(u=0.7,d=0.9,s=2.1),μ=(u=0.3,d=0.2,s=0.6))
        tp=(T=0.8,Φ=0.37,Φbar=0.63,ξ=xi)
        state=(quark_params=qp,thermo_params=tp,K_coeffs=(unused=0.0,))
        pg,pw=AS._cached_interval_grid(0.0,3.4,8)
        extras=custom ? (p_grid=pg,p_w=pw,density_p_grid=pg,density_p_w=pw,
            cos_grid=[-0.8,-0.31,0.12,0.49,0.93],cos_w=[0.1,0.3,0.6,0.7,0.3],
            phi_grid=[0.0,0.6,2pi-0.6,pi,5.1],phi_w=[0.2,0.8,1.1,2.0,0.3]) : (;)
        for process in Main.RelaxTime.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU
            query=s->0.3+0.04s+2exp(-((s-7.1)/0.4)^2)
            plan=prepare_thermal_plan(process,state;p_nodes=8,angle_nodes=5,phi_nodes=9,extras...)
            bounds=ThermalAverageCandidates.process_bounds(process,state;sigma_cutoff=3.05)
            measure=ts.thermal_measure(plan,bounds)
            independent=only(integrate_queries(plan,(query,),(bounds,),Val(:legacy)))
            @test isapprox(ts.integrate_measure(measure,query),independent;rtol=2e-12,atol=1e-14)
        end
    end
end

@testset "Angular geometry is exact and rejects inconsistent inputs" begin
    geometry=reshape([0.0,-0.0,0.1,nextfloat(0.1)],4,1,1)
    @test length(only(AS.angular_groups(geometry,ones(4))).cosines)==4
    @test_throws ArgumentError AS.angular_groups(geometry,ones(3))
    @test_throws ArgumentError AS.angular_groups(fill(NaN,2,1,1),ones(2))
end
