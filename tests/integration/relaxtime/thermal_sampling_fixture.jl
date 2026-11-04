if !isdefined(Main,:RelaxTime)
    Base.include(Main,joinpath(@__DIR__,"..","..","..","src","relaxtime","RelaxTime.jl"))
end

function thermal_sampling_fixture(xi=0.0)
    rt=Main.RelaxTime
    tp=(T=0.9,Φ=0.4,Φbar=0.55,ξ=Float64(xi))
    qp0=(m=(u=0.73,d=0.73,s=1.51),μ=(u=0.2,d=0.2,s=0.5))
    av=rt.AFieldBuilder.build_A_triplet(qp0,tp;cos_nodes=4)
    qp=merge(qp0,(A=av,))
    gu=rt.EffectiveCouplings.calculate_G_from_A(av.u,qp.m.u)
    gs=rt.EffectiveCouplings.calculate_G_from_A(av.s,qp.m.s)
    kc=rt.EffectiveCouplings.calculate_effective_couplings(Main.Constants_PNJL.G_fm2,
        Main.Constants_PNJL.K_fm5,gu,gs)
    options=(p_nodes=4,angle_nodes=3,phi_nodes=5,n_sigma_points=4,
        sampling_config=rt.AverageScatteringRate.ThermalSamplingConfig(seed_order=8))
    return (qp=qp,tp=tp,kc=kc,options=options)
end
