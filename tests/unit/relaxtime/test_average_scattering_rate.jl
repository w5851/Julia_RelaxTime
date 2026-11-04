using Test
if !isdefined(Main,:RelaxTime)
    include(joinpath(@__DIR__,"..","..","..","src","relaxtime","RelaxTime.jl"))
end
if !isdefined(Main,:ThermalAverageCandidates)
    include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","ThermalAverageCandidates.jl"))
end

@testset "Default thermal cache and parameter adapters" begin
    as=Main.RelaxTime.AverageScatteringRate
    rt=Main.RelaxTime
    for xi in (-0.4,0.0,0.4)
        tp=(T=0.9,Φ=0.4,Φbar=0.55,ξ=xi)
        qp0=(m=(u=0.73,d=0.73,s=1.51),μ=(u=0.2,d=0.2,s=0.5))
        av=rt.AFieldBuilder.build_A_triplet(qp0,tp;cos_nodes=4)
        qp=merge(qp0,(A=av,))
        gu=rt.EffectiveCouplings.calculate_G_from_A(av.u,qp.m.u)
        gs=rt.EffectiveCouplings.calculate_G_from_A(av.s,qp.m.s)
        kc=rt.EffectiveCouplings.calculate_effective_couplings(Main.Constants_PNJL.G_fm2,Main.Constants_PNJL.K_fm5,gu,gs)
        options=(p_nodes=4,angle_nodes=3,phi_nodes=5,n_sigma_points=4,sigma_cutoff=2.0,
            sampling_config=as.ThermalSamplingConfig(seed_order=8))
        cache=as.build_thermal_weighted_cache(:us_to_us,qp,tp,kc;options...)
        value=as.average_scattering_rate(:us_to_us,qp,tp,kc;options...,cs_cache=cache)
        qstruct=Main.ParameterTypes.QuarkParams(qp)
        tstruct=Main.ParameterTypes.ThermoParams(tp)
        @test value==cache.rate
        @test as.average_scattering_rate(:us_to_us,qstruct,tstruct,kc;options...,cs_cache=cache)==value
        direct=s->rt.TotalCrossSection.total_cross_section(:us_to_us,s,qp,tp,kc;n_points=4)
        reference=Main.ThermalAverageCandidates.prepared_rate(direct,:us_to_us,
            (quark_params=qp,thermo_params=tp,K_coeffs=kc);
            p_nodes=4,angle_nodes=3,phi_nodes=5,sigma_cutoff=2.0,velocity=:legacy)
        @test isapprox(value,reference;rtol=1e-3)
        for policy in (:default,:validated_anchored)
            @test_throws ArgumentError as.average_scattering_rate(:us_to_us,qp,tp,kc;
                options...,cs_cache=cache,sigma_cache_policy=policy)
        end
    end
end

@testset "Isotropic density and exact RS variable change" begin
    as=Main.RelaxTime.AverageScatteringRate
    p,w=as._cached_semi_infinite_momentum_grid(96,1.0)
    c,cw=as._cached_interval_grid(0.0,1.0,36)
    for species in (:u,:ubar), mu in (0.0,0.4)
        base=as.number_density(species,0.73,mu,0.9,0.4,0.55,0.0;p_grid=p,p_w=w,angle_nodes=36)
        explicit=as.number_density(species,0.73,mu,0.9,0.4,0.55,0.0;
            p_grid=p,p_w=w,cos_grid=c,cos_w=cw,angle_nodes=36)
        @test base==explicit
        for xi in (-0.5,0.0,0.5)
            density=as.number_density(species,0.73,mu,0.9,0.4,0.55,xi;p_grid=p,p_w=w,angle_nodes=36)
            @test isapprox(density,base/sqrt(1+xi);rtol=1e-7)
        end
    end
end
