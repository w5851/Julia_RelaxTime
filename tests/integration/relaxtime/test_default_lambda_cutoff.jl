using Test
include(joinpath(@__DIR__,"thermal_sampling_fixture.jl"))

@testset "Thermal cache keeps the physical Lambda cutoff" begin
    as=Main.RelaxTime.AverageScatteringRate
    rt=Main.RelaxTime.RelaxationTime
    f=thermal_sampling_fixture(0.2)
    cutoff=Main.Constants_PNJL.Λ_inv_fm
    process=:ssbar_to_uubar
    cache=as.build_thermal_weighted_cache(process,f.qp,f.tp,f.kc;f.options...)
    explicit=as.build_thermal_weighted_cache(process,f.qp,f.tp,f.kc;f.options...,sigma_cutoff=cutoff)
    @test cache.identity==explicit.identity
    @test cache.rate==explicit.rate
    @test last(cache.context.bounds)==(2sqrt(f.qp.m.u^2+cutoff^2))^2
    # Reuse the heavy-initial channel through the full tau composition path.
    existing=Dict(p=>0.5 for p in rt.REQUIRED_RATE_KEYS_FOR_TAU if p!=process)
    densities=(u=.6,d=.6,s=.2,ubar=.3,dbar=.3,sbar=.1)
    options=(f.options...,densities=densities,existing_rates=existing,cs_caches=Dict(process=>cache))
    automatic=rt.relaxation_times(f.qp,f.tp,f.kc;options...)
    manual=rt.relaxation_times(f.qp,f.tp,f.kc;options...,sigma_cutoff=cutoff)
    @test automatic.tau==manual.tau
    @test automatic.rates.ssbar_to_uubar==cache.rate
end
