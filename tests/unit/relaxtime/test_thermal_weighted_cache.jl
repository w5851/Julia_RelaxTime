using Test
if !isdefined(Main,:RelaxTime)
    include(joinpath(@__DIR__,"..","..","..","src","relaxtime","RelaxTime.jl"))
end
if !isdefined(Main,:ThermalAverageCandidates)
    include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","ThermalAverageCandidates.jl"))
end

@testset "Thermal cache binds physical inputs, quadrature and accuracy" begin
    as=Main.RelaxTime.AverageScatteringRate
    rt=Main.RelaxTime
    tp=(T=0.9,Φ=0.4,Φbar=0.55,ξ=0.27)
    qp0=(m=(u=0.73,d=0.73,s=1.51),μ=(u=0.2,d=0.2,s=0.5))
    av=rt.AFieldBuilder.build_A_triplet(qp0,tp;cos_nodes=4)
    qp=merge(qp0,(A=av,))
    gu=rt.EffectiveCouplings.calculate_G_from_A(av.u,qp.m.u)
    gs=rt.EffectiveCouplings.calculate_G_from_A(av.s,qp.m.s)
    kc=rt.EffectiveCouplings.calculate_effective_couplings(Main.Constants_PNJL.G_fm2,
        Main.Constants_PNJL.K_fm5,gu,gs)
    process=:us_to_us
    cfg=as.ThermalSamplingConfig(rtol=1e-3,seed_order=8)
    options=(p_nodes=4,angle_nodes=3,phi_nodes=5,n_sigma_points=4,sigma_cutoff=2.0,sampling_config=cfg)
    cache=as.build_thermal_weighted_cache(process,qp,tp,kc;options...)
    @test cache.statistics.guard_phase==2
    @test cache.estimated_error<=cfg.atol+cfg.rtol*cache.rate
    @test cache.rate>0
    @test cache.identity.sampling.relative_velocity === :cm
    @test as.ThermalSamplingConfig(relative_velocity="medium_moller").relative_velocity === :medium_moller
    @test as.average_scattering_rate(process,qp,tp,kc;options...,cs_cache=cache,
        sigma_cache_policy=:thermal_weighted)==cache.rate
    # The outer A field is not used by the distributions; the propagator A is
    # rebuilt with the same documented defaults when the caller omits it.
    @test as.average_scattering_rate(process,qp0,tp,kc;options...,cs_cache=cache,
        sigma_cache_policy=:thermal_weighted)==cache.rate

    state=(quark_params=qp,thermo_params=tp,K_coeffs=kc)
    direct=s->rt.TotalCrossSection.total_cross_section(process,s,qp,tp,kc;n_points=4)
    reference=Main.ThermalAverageCandidates.prepared_rate(direct,process,state;
        p_nodes=4,angle_nodes=3,phi_nodes=5,sigma_cutoff=2.0,velocity=:legacy)
    @test isapprox(cache.rate,reference;rtol=1e-3)
    middle=cache.table.nodes[length(cache.table.nodes)÷2]
    @test as.get_sigma(cache,middle)==cache.table(middle)
    @test as.get_sigma(cache,middle;direct=true)==direct(middle)
    @test as.get_sigma(cache,first(cache.context.bounds))==0
    @test as.get_sigma(cache,last(cache.context.bounds))==0
    before=copy(cache.table.sigma)
    outside=(first(cache.context.bounds)+first(cache.table.nodes))/2
    @test_throws DomainError as.get_sigma(cache,outside)
    @test as.get_sigma(cache,outside;direct=true)>=0
    @test cache.table.sigma==before

    for changed in ((p_nodes=5,),(angle_nodes=4,),(phi_nodes=6,),
        (density_p_nodes=48,),(scale=5.0,),(density_scale=5.0,),
        (n_sigma_points=6,),(sigma_cutoff=2.1,),
        (sampling_config=as.ThermalSamplingConfig(rtol=5e-4,seed_order=8),))
        @test_throws ArgumentError as.average_scattering_rate(process,qp,tp,kc;
            merge(options,changed)...,cs_cache=cache,sigma_cache_policy=:thermal_weighted)
    end
    for changed in (merge(tp,(T=0.91,)),merge(tp,(ξ=0.28,)),merge(tp,(Φbar=0.54,)))
        @test_throws ArgumentError as.average_scattering_rate(process,qp,changed,kc;
            options...,cs_cache=cache,sigma_cache_policy=:thermal_weighted)
    end
    changed_q=merge(qp,(A=merge(qp.A,(u=qp.A.u+1e-3,)),))
    @test_throws ArgumentError as.average_scattering_rate(process,changed_q,tp,kc;
        options...,cs_cache=cache,sigma_cache_policy=:thermal_weighted)
    @test_throws ArgumentError as.average_scattering_rate(:ud_to_ud,qp,tp,kc;
        options...,cs_cache=cache,sigma_cache_policy=:thermal_weighted)
    @test as.average_scattering_rate(process,qp,tp,kc;options...,cs_cache=cache)==cache.rate
    @test_throws ArgumentError as.average_scattering_rate(process,qp,tp,kc;
        merge(options,(sampling_config=as.ThermalSamplingConfig(seed_order=8,
            relative_velocity=:medium_moller),))...,
        cs_cache=cache,sigma_cache_policy=:thermal_weighted)
    @test_throws ArgumentError as.average_scattering_rate(process,qp,tp,kc;
        options...,cs_cache=cache,sigma_cache_policy=:validated_anchored)
    @test_throws ArgumentError as.average_scattering_rate(process,qp,tp,kc;
        options...,cs_cache=cache,sigma_cache_policy=:thermal_weighted,interpolation_mode=:linear)
    @test_throws ArgumentError as.average_scattering_rate(process,qp,tp,kc;
        options...,cs_cache=cache,sigma_cache_policy=:thermal_weighted,angular_reduction=:none)
    @test_throws ArgumentError as.build_thermal_weighted_cache(process,qp,tp,kc;
        merge(options,(sampling_config=as.ThermalSamplingConfig(max_evaluations=1),))...)
    @test_throws ArgumentError as.build_thermal_weighted_cache(process,qp,tp,kc;
        merge(options,(sampling_config=as.ThermalSamplingConfig(max_events=1),))...)
    @test_throws ArgumentError as.ThermalSamplingConfig(rtol=0,atol=0)
    @test_throws ArgumentError as.ThermalSamplingConfig(max_rounds=0)

    # Reuse through the actual relaxation-time entrypoint, with one missing
    # channel to ensure the new cache is consumed rather than bypassed.
    existing=Dict(p=>0.5 for p in rt.RelaxationTime.REQUIRED_RATE_KEYS_FOR_TAU if p!=process)
    densities=(u=0.6,d=0.6,s=0.2,ubar=0.3,dbar=0.3,sbar=0.1)
    answer=rt.RelaxationTime.relaxation_times(qp,tp,kc;densities=densities,
        existing_rates=existing,cs_caches=Dict(process=>cache),options...,
        sigma_cache_policy=:thermal_weighted)
    @test answer.rates.us_to_us==cache.rate
    @test answer.tau.u==1/rt.RelaxationTime.relaxation_rates(densities,answer.rates).u
end
